# import copy
import os
from pdb_processing import *
from Bio.PDB import Superimposer, PDBIO



def testfunc(out_dir):
    # Get the ModelChain with the longest interactions list
    first_modelchain = max(processed_chains, key = lambda x: len(x.interactions))

    # Objects with the first chain to use and the first Structure (all its chains) to use
    first_chain = first_modelchain.interactions[0][0]
    first_pdb = first_modelchain.interactions[0][0].get_parent().copy()

    # The first interactor: the second element in the interactions list
    # (the first element's interaction is already in the full Structure first_pdb)
    first_interaction = first_modelchain.interactions[1]

    # Superimpose the whole Structure of the first interaction into the first_pdb
    # The superimposition has to be between the specific chains
    ref_atoms = list(first_chain.get_atoms())

    # The movable chain that belongs to the same ModelChain as first_chain
    mov_atoms = list(first_interaction[0].get_atoms())

    # Get the same amount of atoms
    if len(ref_atoms) != len(mov_atoms):
        atoms_length = min(len(ref_atoms), len(mov_atoms))
        mov_atoms_idx = (len(mov_atoms) - atoms_length)/2
        mov_atoms = mov_atoms[mov_atoms_idx:len(mov_atoms)-mov_atoms_idx]
        ref_atoms_idx = (len(ref_atoms) - atoms_length)/2
        ref_atoms = ref_atoms[ref_atoms_idx:len(ref_atoms)-ref_atoms_idx]

    # Perform the superimposition and save the translation and rotation vectors
    superimpos = Superimposer()
    superimpos.set_atoms(ref_atoms, mov_atoms)
    mov = superimpos.rotran

    # Make a copy of the chain that has to move:
    # From the second interactor tuple, the last element (which is
    # the one the superimposed chain interacts with)
    first_interactor = first_interaction[-1].copy()
    first_interactor.transform(mov[0], mov[1])
    first_pdb.add(first_interactor)

    # Save the PDB of the complex
    io = PDBIO()
    io.set_structure(first_pdb)
    pdb_name = "test"
    io.save(os.path.join(out_dir, "structures", pdb_name + ".pdb"))
