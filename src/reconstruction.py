# import copy
import os
from pdb_processing import *
from Bio.PDB import Superimposer, PDBIO
from string import ascii_uppercase

letter_list = list(ascii_uppercase)  # List of A-Z to use as the chains' new ids
chain_ids = letter_list + [a + b for a in letter_list for b in letter_list]  # Extended list A-ZZ

added_chains = dict()


def testfunc(out_dir):
    # Get the ModelChain with the longest interactions list
    first_modelchain = max(processed_chains, key = lambda x: len(x.interactions))

    # Objects with the first chain to use and the first Structure (all its chains) to use
    first_chain = first_modelchain.interactions[0][0]
    first_pdb = first_modelchain.interactions[0][0].parent.parent.copy()

    for chain in list(first_pdb.get_chains()):
        added_chain_id(chain)

    # "Mark" the chain as processed so that its interactions aren't tried to be added a second time
    list(first_pdb.get_chains())[0].xtra = "processed"
    add_modelchain_interactions(first_pdb, first_chain, first_modelchain)



    for chain in [chain for chain in list(first_pdb.get_chains()) if chain.xtra != "processed"]:
        chain_full_id = added_chains[chain.get_id()]
        modelchain_obj = chain_to_model_chain[chain_full_id]


    # Save the PDB of the complex
    save_pdb(first_pdb, out_dir)




def added_chain_id(chain):
    new_id = chain_ids.pop(0)
    chain_full_id = get_chain_full_id(chain)
    added_chains[new_id] = chain_full_id
    chain.id = new_id
    chain.xtra = chain_full_id
    return 0



def add_modelchain_interactions(structure, ref_chain, modelchain_obj):
    for interaction in modelchain_obj.interactions:
        mov = get_rotran_matrix(ref_chain, interaction[0])
        interactor = interaction[-1].copy()
        added_chain_id(interactor)
        interactor.transform(mov[0], mov[1])
        if not is_clashing(structure, interactor):
            structure[0].add(interactor)
        else:
            continue


def get_rotran_matrix(ref_chain, mov_chain):
    # Superimpose the whole Structure of the first interaction into the first_pdb
    # The superimposition has to be between the specific chains
    ref_atoms = list(ref_chain.get_atoms())

    # The movable chain that belongs to the same ModelChain as first_chain
    mov_atoms = list(mov_chain.get_atoms())

    # Get the same amount of atoms
    if len(ref_atoms) != len(mov_atoms):
        atoms_length = min(len(ref_atoms), len(mov_atoms))
        ref_atoms_idx, remainder = divmod((len(ref_atoms) - atoms_length),  2)
        ref_atoms = ref_atoms[ref_atoms_idx:len(ref_atoms) - ref_atoms_idx - remainder]
        mov_atoms_idx, remainder = divmod((len(mov_atoms) - atoms_length),  2)
        mov_atoms = mov_atoms[mov_atoms_idx:len(mov_atoms) - mov_atoms_idx - remainder]

    # Perform the superimposition and save the translation and rotation vectors
    superimpos = Superimposer()
    superimpos.set_atoms(ref_atoms, mov_atoms)
    return superimpos.rotran


def is_clashing(structure, interactor):
    atoms_int = list(interactor.get_atoms())
    atoms_str = list(structure.get_atoms())

    neighbors_search = NeighborSearch(atoms_int)
    clashes = 0
    for atom in atoms_str:
        close_atoms = neighbors_search.search(atom.coord, 0.5)
        clashes += len(close_atoms)
    if clashes < 10:
        return False
    else:
        return True







def save_pdb(structure, out_dir):
    io = PDBIO()
    io.set_structure(structure)
    pdb_name = "Complex"
    io.save(os.path.join(out_dir, "structures", pdb_name + ".pdb"))
