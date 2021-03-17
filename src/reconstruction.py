# import copy
import os
from pdb_processing import *
from Bio.PDB import Superimposer, PDBIO
from string import ascii_uppercase

letter_list = list(ascii_uppercase)  # List of A-Z to use as the chains' new ids
chain_ids = letter_list + [a + b for a in letter_list for b in letter_list]  # Extended list A-ZZ

added_chains = dict()


def build_complex(out_dir):
    # Get the ModelChain with the longest interactions list
    first_modelchain = max(processed_chains, key = lambda x: len(x.interactions))

    # Initialize the complex Structure object with one of the PDBs of the ModelChain
    first_chain = first_modelchain.interactions[0][0]
    complex = first_chain.parent.parent.copy()

    # Process the first chains: add them to the added_chains dictionary and change their ID
    for chain in list(complex.get_chains()):
        rename_added_chain(chain)

    chain_number_change = 1
    # Keep trying to add new interactions until no new chains are added
    while chain_number_change > 0:
        # Make a list of the chains in the complex and save the number of chains
        chain_list = list(complex.get_chains())
        chain_number = len(chain_list)
        print("Chains before adding:", chain_number)

        # For each chain in the complex whose interactions haven't already been added (not
        # "processed"), add all the interactions saved in its corresponding ModelChain
        for chain in [chain for chain in chain_list if "processed" not in chain.xtra]:
            # Retrieve the corresponding ModelChain by querying the chain_full_id in the
            # added_chains dictionary and the ModelChain object in the chain_to_model_chain dict
            chain_full_id = added_chains[chain.get_id()]
            modelchain_obj = chain_to_model_chain[chain_full_id]
            print("Adding interactions of:", chain_full_id, len(modelchain_obj.interactions))

            # Add all the chains that interact that are present in the ModelChain, if they
            # don't clash with existing chains of the complex
            add_modelchain_interactions(complex, chain, modelchain_obj)
            # Change the xtra attribute of the chain whose interactions have been added so that
            # they aren't tried to be added again
            chain.xtra["processed"] = True

        # Get the new number of chains in the complex to keep adding interactors if new chains have
        # been added or stop if none have been added or we have reached the max number of chains
        new_chain_number = len(list(complex.get_chains()))
        print("Chains after adding:", new_chain_number)
        if new_chain_number <= 100:
            chain_number_change = new_chain_number - chain_number
        else:
            chain_number_change = 0


    # Finally, save the PDB of the Complex
    save_pdb(complex, out_dir)
    return 0

def save_pdb(structure, out_dir):
    io = PDBIO()
    io.set_structure(structure)
    pdb_name = "Complex"
    io.save(os.path.join(out_dir, "structures", pdb_name + ".pdb"))



def rename_added_chain(chain):
    # Get the new ID for the chain from the chain_ids list
    new_id = chain_ids.pop(0)
    # Get the full ID of the chain to map it to its new chain ID with the added_chains dictionary
    chain_full_id = get_chain_full_id(chain)
    # Save the new information
    added_chains[new_id] = chain_full_id
    chain.id = new_id
    chain.xtra["full_id"] = chain_full_id
    return 0



def add_modelchain_interactions(structure, ref_chain, modelchain_obj):
    # For every interaction in the ModelChain
    for interaction in modelchain_obj.interactions:
        # Get the rotation-translation matrix
        mov = get_rotran_matrix(ref_chain, interaction[0])
        # Make a copy of the chain to add it to the Complex being built
        interactor = interaction[-1].copy()
        # Apply the rotation-translation
        interactor.transform(mov[0], mov[1])
        # If its new atom coordinates won't clash with any existing atoms in the Complex, add it
        if not is_clashing(structure, interactor):
            # Process the chain IDs before adding it
            rename_added_chain(interactor)
            structure[0].add(interactor)
        else:
            continue
    return 0


def get_rotran_matrix(ref_chain, mov_chain):
    ref_atoms = list(ref_chain.get_atoms())
    mov_atoms = list(mov_chain.get_atoms())

    # Get the same amount of atoms
    if len(ref_atoms) != len(mov_atoms):
        atoms_length = min(len(ref_atoms), len(mov_atoms))
        ref_atoms_idx, remainder = divmod((len(ref_atoms) - atoms_length),  2)
        ref_atoms = ref_atoms[ref_atoms_idx:len(ref_atoms) - ref_atoms_idx - remainder]
        mov_atoms_idx, remainder = divmod((len(mov_atoms) - atoms_length),  2)
        mov_atoms = mov_atoms[mov_atoms_idx:len(mov_atoms) - mov_atoms_idx - remainder]

    # Perform the superimposition and save the translation-rotation matrix
    superimpos = Superimposer()
    superimpos.set_atoms(ref_atoms, mov_atoms)
    return superimpos.rotran


def is_clashing(structure, interactor):
    atoms_int = list(interactor.get_atoms())
    atoms_str = list(structure.get_atoms())

    neighbors_search = NeighborSearch(atoms_int)
    clashes = 0
    # For each atom in the Complex, check if there are clashing atoms in the new chain
    for atom in atoms_str:
        close_atoms = neighbors_search.search(atom.coord, 0.01)
        clashes += len(close_atoms)
    # Don't consider clashing if there are less than a certain number of close atoms
    if clashes < 10:
        return False
    else:
        return True
