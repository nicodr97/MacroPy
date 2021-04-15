import logging as log
from PDB_processing import processed_chains, chain_to_model_chain
from Bio.PDB import NeighborSearch, Structure, Model
from string import ascii_uppercase
from PDB_tools import align, get_chain_full_id, save_structure, superimpose
from Stoichiometry import get_common_chain_id, is_in_stoichiometry, update_stoichiometry

letter_list = list(ascii_uppercase)  # List of A-Z to use as the chains' new ids
chain_ids = letter_list + [a + b for a in letter_list for b in letter_list]  # Extended list A-ZZ
current_stoich_dict = dict()  # Count the stoichiometry in the reconstructed complex
current_stoich_dict_by_prefix = dict()  # Count the chain stoichiometry when the stoichiometry


# file has Uniprot prefixes as ids


def build_complex(out_dir, max_chains, number_clashes, save_pdb, clashes_distance,
                  stoich_dict, complex_name, identity_threshold):
    """Main function of the process of building the complex"""

    log.info("Building new complex")
    first_modelchain = choose_first_modelchain(stoich_dict)

    # Initialize the complex Structure object with one of the PDBs of the ModelChain
    first_chain = first_modelchain.chain
    macro_complex = Structure.Structure(complex_name)
    macro_complex.add(Model.Model(0))
    macro_complex[0].add(first_chain.copy())
    # Save the original ID and change its ID in the complex
    macro_complex[0][first_chain.get_id()].xtra["full_id"] = get_chain_full_id(first_chain)
    macro_complex[0][first_chain.get_id()].id = chain_ids.pop(0)
    log.info(f"First chain in the complex: {get_chain_full_id(first_chain)}")

    if stoich_dict:
        update_stoichiometry(stoich_dict, first_chain, current_stoich_dict,
                             current_stoich_dict_by_prefix)

    chain_number_change = 1
    # Keep trying to add new interactions until no new chains are added
    while chain_number_change > 0 and len(list(macro_complex.get_chains())) <= max_chains:
        # Make a list of the chains in the complex and save the number of chains
        chain_list = list(macro_complex.get_chains())
        chain_number = len(chain_list)
        log.info(f"Chains in the complex: {chain_number}")

        # For each chain in the complex whose interactions haven't already been added (not
        # "processed"), add all the interactions saved in its corresponding ModelChain
        for chain in [chain for chain in chain_list if "processed" not in chain.xtra]:
            # Retrieve the corresponding ModelChain by obtaining the chain_full_id from the .xtra
            # attribute and querying the chain_to_model_chain dict for the ModelChain object
            chain_full_id = chain.xtra["full_id"]
            modelchain_obj = chain_to_model_chain[chain_full_id]
            log.info(f"Checking the {len(modelchain_obj.interactions)} interactions of chain"
                     f" {chain_full_id}")

            # Add all the chains that interact that are present in the ModelChain, if they
            # don't clash with existing chains of the complex
            add_modelchain_interactions(macro_complex, chain, modelchain_obj, clashes_distance,
                                        number_clashes, stoich_dict, identity_threshold)
            # Change the xtra attribute of the chain whose interactions have been added so that
            # they aren't tried to be added again
            chain.xtra["processed"] = True

        # Get the new number of chains in the complex to keep adding interactors if new chains have
        # been added or stop if none have been added or we have reached the max number of chains
        new_chain_number = len(list(macro_complex.get_chains()))
        chain_number_change = new_chain_number - chain_number

    total_chains = len(list(macro_complex.get_chains()))
    log.info(f"The final complex has {total_chains} chains")
    if current_stoich_dict != stoich_dict:
        log.error("Stoichiometry couldn't be satisfied")
    # Finally, save the PDB of the Complex
    save_structure(macro_complex, out_dir, complex_name, save_pdb)


def choose_first_modelchain(stoich_dict):
    """
    Select starting chain of the complex based on:
        - The nucleotide chain with longest sequence
        - The protein chain with more interactions otherwise
    """

    if any(model_chain.chain.xtra["type"] == "nuc" and (
            not stoich_dict or is_in_stoichiometry(stoich_dict, model_chain.chain)) for model_chain
           in processed_chains):
        return max(filter(lambda x: x.chain.xtra["type"] == "nuc" and (
                not stoich_dict or is_in_stoichiometry(stoich_dict, x.chain)),
                          processed_chains), key=lambda x: len(x.sequence))

    return max([model_chain for model_chain in processed_chains
                if not stoich_dict or is_in_stoichiometry(stoich_dict, model_chain.chain)],
               key=lambda x: len(set(get_chain_full_id(inter[0]).split("_")[0]
                                     for inter in x.interactions)))


def add_modelchain_interactions(structure, ref_chain, modelchain_obj, clashes_distance,
                                number_clashes, stoich_dict, identity_threshold):
    """Add the chains interacting with the reference chain based on certain restrictions"""

    for interaction in ordered_interactions(modelchain_obj, ref_chain):
        # Continue if there isn't a stoichiometry file or, if there is, if the chain is in it
        if not stoich_dict or is_in_stoichiometry(stoich_dict, interaction[-1]):
            # Make a copy of the chain to add it to the Complex being built
            interactor = interaction[-1].copy()

            if stoich_dict:
                common_chain_id = get_common_chain_id(stoich_dict, interactor)
                if common_chain_id in current_stoich_dict and current_stoich_dict[
                    common_chain_id] >= stoich_dict[common_chain_id]:
                    # Skip interaction if the stoichiometry is exceeded
                    continue

            # Add the interaction of the ref_chain to every position it aligns to (e.g., for a
            # long nucleic acid chain, every time a shorter nucleic acid chain aligns with it)
            for percentage, positions in align(ref_chain, interaction[0].xtra["seq"]):
                # Add the interaction if the alignment score % is higher than the threshold
                if percentage > identity_threshold:
                    # Get the rotation-translation matrix
                    rmsd, mov = superimpose(ref_chain, interaction[0], positions)
                    # Apply the rotation-translation
                    interactor.transform(mov[0], mov[1])
                    if not is_clashing(structure, interactor, clashes_distance, number_clashes):
                        if stoich_dict:
                            update_stoichiometry(stoich_dict, interaction[-1], current_stoich_dict,
                                                 current_stoich_dict_by_prefix)
                        # Process the chain IDs before adding it
                        interactor.xtra["full_id"] = get_chain_full_id(interactor)
                        interactor.id = chain_ids.pop(0)
                        structure[0].add(interactor)
                        log.info(f"Added chain {interactor.xtra['full_id']}")


def ordered_interactions(modelchain_obj, ref_chain):
    """Return a generator of the interactions of a chain sorted by length
    if the chain is a nucleotide and has a complementary strand"""

    # If the chain to which we are adding interactions is a nucleic acid
    # and at least 1 (any) of the chains it interacts with is also a nucleic acid
    if ref_chain.xtra["type"] == "nuc" and any([int[-1].xtra["type"] == "nuc"
                                                for int in modelchain_obj.interactions]):
        # Get the name of the PDB file (e.g., 2v2t) of the chain
        chain_pdb_name = ref_chain.xtra["full_id"].split(".")[-1][:4]
        # Get the list of interactions with a nucleic acid chain that share the same PDB file
        chains_list = [int for int in modelchain_obj.interactions
                       if int[-1].xtra["type"] == "nuc" and
                       int[0].xtra["full_id"].split(".")[-1][:4] == chain_pdb_name]

        if len(chains_list) > 0:
            # Get the index of the interactions list to start iterating from
            index = modelchain_obj.interactions.index(chains_list[0])
            for i in modelchain_obj.interactions[index:]:
                yield i
            if index != 0:
                for i in modelchain_obj.interactions[:index]:
                    yield i
    # If the chain to which we are adding interactions is a protein, proceed normally
    else:
        for i in modelchain_obj.interactions:
            yield i


def is_clashing(structure, interactor, clashes_distance, number_clashes):
    """
    Check if two chains (structure and interactor) are clashing,
    i.e. if they have more than a threshold of clashes
    """

    # First, search for close CAs and Ps to then restrict an exhaustive search to the closer chains
    ca_str = [atom for atom in structure.get_atoms() if
              atom.get_id() == "CA" or atom.get_id() == "P"]
    ca_int = [atom for atom in interactor.get_atoms() if
              atom.get_id() == "CA" or atom.get_id() == "P"]

    neighbors_search = NeighborSearch(ca_str)
    close_ca = list()
    # For each CA or P in the interactor, check if there are close CAs or Ps in the complex
    for atom in ca_int:
        close_atoms = neighbors_search.search(atom.coord, 5)
        if len(close_atoms) > 0:
            close_ca.append(close_atoms)

    # If there are none, don't consider clashing
    if len(close_ca) == 0:
        return False

    # Get the chains of the Complex from which at least a CA has been found to be close
    close_chains = set(atom.parent.parent for ca_list in close_ca for atom in ca_list)

    # Exhaustive search, with all the atoms from the closer chains and the interactor
    chains_atoms = [atom for chain in close_chains for atom in chain.get_atoms()]
    int_atoms = list(interactor.get_atoms())

    exhaustive_search = NeighborSearch(chains_atoms)
    clashes = 0
    # For each atom in the interactor, check if there are clashing atoms in the Complex
    for atom in int_atoms:
        close_atoms = exhaustive_search.search(atom.coord, clashes_distance)
        clashes += len(close_atoms)
        if clashes >= number_clashes:
            return True
    # Don't consider clashing if there are less than a certain number of close atoms
    return False
