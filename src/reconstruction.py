import os
import sys
import logging as log
from pdb_processing import processed_chains, chain_to_model_chain, get_chain_full_id, \
    get_atom_chains_with_same_length, model_chain_to_chains
from Bio.PDB import Superimposer, MMCIFIO, NeighborSearch, Structure, Model
from string import ascii_letters

chain_ids = list(ascii_letters) + [str(x) for x in range(150)]  # List of A-Z, a-z and numbers to
# use as the chains' new ids
MAX_CHAINS = 180
current_stoich_dict = dict()  # Count the stoichiometry of each chain in the complex


def build_complex(out_dir, clashes_distance, ca_distance, number_clashes, stoich_dict):
    # Get the ModelChain with the longest interactions list
    # first_modelchain = max(processed_chains, key = lambda x: len(x.interactions))
    first_modelchain = choose_first_modelchain(stoich_dict)

    # Initialize the complex Structure object with one of the PDBs of the ModelChain
    first_chain = first_modelchain.chain
    macro_complex = Structure.Structure(first_chain.parent.parent.get_id())
    macro_complex.add(Model.Model(0))
    macro_complex[0].add(first_chain.copy())
    if stoich_dict:
        current_stoich_dict[get_common_chain_id(stoich_dict, first_chain)] = 1

    # Process the first chain: change its ID and save the original full ID
    rename_added_chain(macro_complex[0][first_chain.get_id()])

    chain_number_change = 1
    # Keep trying to add new interactions until no new chains are added
    while chain_number_change > 0 and len(list(macro_complex.get_chains())) <= MAX_CHAINS:
        # Make a list of the chains in the complex and save the number of chains
        chain_list = list(macro_complex.get_chains())
        chain_number = len(chain_list)
        log.info(f"Chains before adding: {chain_number}")

        # For each chain in the complex whose interactions haven't already been added (not
        # "processed"), add all the interactions saved in its corresponding ModelChain
        for chain in [chain for chain in chain_list if "processed" not in chain.xtra]:
            # Retrieve the corresponding ModelChain by obtaining the chain_full_id from the .xtra
            # attribute and querying the chain_to_model_chain dict for the ModelChain object
            chain_full_id = chain.xtra["full_id"]
            modelchain_obj = chain_to_model_chain[chain_full_id]
            log.info(f"Adding interactions of: {chain_full_id} {len(modelchain_obj.interactions)}")

            # Add all the chains that interact that are present in the ModelChain, if they
            # don't clash with existing chains of the complex
            add_modelchain_interactions(macro_complex, chain, modelchain_obj, clashes_distance,
                                        ca_distance, number_clashes, stoich_dict)
            # Change the xtra attribute of the chain whose interactions have been added so that
            # they aren't tried to be added again
            chain.xtra["processed"] = True

        # Get the new number of chains in the complex to keep adding interactors if new chains have
        # been added or stop if none have been added or we have reached the max number of chains
        new_chain_number = len(list(macro_complex.get_chains()))
        log.info(f"Chains after adding: {new_chain_number}")
        chain_number_change = new_chain_number - chain_number

    # Finally, save the PDB of the Complex
    save_pdb(macro_complex, out_dir)


def save_pdb(structure, out_dir):
    io = MMCIFIO()
    io.set_structure(structure)
    pdb_name = "Complex"
    io.save(os.path.join(out_dir, "structures", pdb_name + ".cif"))


def choose_first_modelchain(stoich_dict):
    # if any(model_chain.chain.xtra["type"] == "nuc" for model_chain in processed_chains):
    #     return max(filter(lambda x: x.chain.xtra["type"] == "nuc", processed_chains), key=lambda x: len(x.sequence))

    return max([model_chain for model_chain in processed_chains
                if not stoich_dict or check_stoichiometry(stoich_dict, model_chain.chain)],
               key=lambda x: len(set(get_chain_full_id(inter[0]).split("_")[0]
                                                       for inter in x.interactions)))


def rename_added_chain(chain):
    # Get the new ID for the chain from the chain_ids list
    new_id = chain_ids.pop(0)
    # Get the full ID of the chain to save it in the .xtra attribute
    chain_full_id = get_chain_full_id(chain)
    # Save the new information
    chain.id = new_id
    chain.xtra["full_id"] = chain_full_id


def add_modelchain_interactions(structure, ref_chain, modelchain_obj, clashes_distance,
                                ca_distance, number_clashes, stoich_dict):
    # For every interaction in the ModelChain
    for interaction in modelchain_obj.interactions:
        if not stoich_dict or check_stoichiometry(stoich_dict, interaction[-1]):
        # Make a copy of the chain to add it to the Complex being built
            interactor = interaction[-1].copy()

            if stoich_dict:
                common_chain_id = get_common_chain_id(stoich_dict, interactor)
                if common_chain_id in current_stoich_dict and current_stoich_dict[
                    common_chain_id] >= stoich_dict[common_chain_id]:
                    # Skip interaction if the stoichiometry is exceeded
                    continue

            # Get the rotation-translation matrix
            mov = get_rotran_matrix(ref_chain, interaction[0])
        # Apply the rotation-translation
            interactor.transform(mov[0], mov[1])
        # If its new atom coordinates won't clash with any existing atoms in the Complex, add it
            if not is_clashing(structure, interactor, clashes_distance, ca_distance,
                               number_clashes):
                if stoich_dict:
                    current_stoich_dict[common_chain_id] = current_stoich_dict.setdefault(
                        common_chain_id, 0) + 1
            # Process the chain IDs before adding it
                rename_added_chain(interactor)
                structure[0].add(interactor)


def get_rotran_matrix(ref_chain, mov_chain):
    ref_atoms = list(ref_chain.get_atoms())
    mov_atoms = list(mov_chain.get_atoms())

    # Get the same amount of atoms
    ref_atoms, mov_atoms = get_atom_chains_with_same_length(ref_atoms, mov_atoms)

    # Perform the superimposition and save the translation-rotation matrix
    superimposer = Superimposer()
    superimposer.set_atoms(ref_atoms, mov_atoms)
    return superimposer.rotran


def is_clashing(structure, interactor, clashes_distance, ca_distance, number_clashes):
    # First, search for close CAs to then restrict a more exhaustive search to the closer chains
    ca_str = [atom for atom in structure.get_atoms() if atom.get_id() == "CA" or
              atom.get_id() == "P"]
    ca_int = [atom for atom in interactor.get_atoms() if atom.get_id() == "CA" or
              atom.get_id() == "P"]

    neighbors_search = NeighborSearch(ca_str)
    close_ca = list()
    # For each CA in the interactor, check if there are close CAs in the complex
    for atom in ca_int:
        close_atoms = neighbors_search.search(atom.coord,
                                              float(clashes_distance) * float(ca_distance))
        if len(close_atoms) > 0:
            close_ca.append(close_atoms)

    # Get the chains of the Complex from which at least a CA has been found to be close
    close_chains = set(atom.parent.parent for ca_list in close_ca for atom in ca_list)
    # If there are none, don't consider clashing
    if len(close_chains) == 0:
        return False

    # Exhaustive search, with all the atoms from the closer chains and the interactor
    chains_atoms = [atom for chain in close_chains for atom in chain.get_atoms()]
    int_atoms = list(interactor.get_atoms())

    exhaustive_search = NeighborSearch(chains_atoms)
    clashes = 0
    # For each atom in the interactor, check if there are clashing atoms in the Complex
    for atom in int_atoms:
        close_atoms = exhaustive_search.search(atom.coord, float(clashes_distance))
        clashes += len(close_atoms)
        if clashes >= int(number_clashes):
            return True
    # Don't consider clashing if there are less than a certain number of close atoms
    return False


# Check if the chain is present in the stoichiometry
def check_stoichiometry(stoich_dict, chain):
    model_chain = chain_to_model_chain[get_chain_full_id(chain)]
    chains = model_chain_to_chains[model_chain.id]
    return any([chain in stoich_dict for chain in chains])


# Get the chain id used in the stoichiometry to use as common chain id
def get_common_chain_id(stoich_dict, chain):
    model_chain = chain_to_model_chain[get_chain_full_id(chain)]
    chains = model_chain_to_chains[model_chain.id]
    common_chain_id = [chain for chain in chains if chain in stoich_dict]
    if len(common_chain_id) > 1:
        log.error(f"Error with stoichiometry: {common_chain_id} are the same chains but different "
                  "lines in the stoichiometry file")
        sys.exit(1)

    return common_chain_id[0]
