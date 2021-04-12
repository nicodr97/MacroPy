import os
import logging as log
from Bio.Align import PairwiseAligner
from Bio.PDB import Superimposer, MMCIFIO



def get_chain_full_id(chain):
    return ":".join(chain.get_full_id()[0::2])




def align(chain, seq):
    # Get chain sequence from the xtra attribute
    chain_seq = chain.xtra["seq"]
    # Compare the chain_seq with the sequence (seq)
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    alignment = aligner.align(chain_seq, seq)

    # Yield all possible alignments found by the function
    for i in alignment:
        lengths = (len(chain_seq), len(seq))
        length = max(lengths) if chain.xtra["type"] == "prot" else min(lengths)
        yield (i.score/length), i.aligned




def superimpose(ref_chain, mov_chain, alignment_positions=None):
    if alignment_positions is not None:
        start_res, end_res = alignment_positions[0][0]
        chain_res = list(ref_chain.get_residues())[start_res:end_res]
        ref_atoms = [atom for res in chain_res for atom in res.get_atoms()]

        start_res, end_res = alignment_positions[1][0]
        mov_res = list(mov_chain.get_residues())[start_res:end_res]
        mov_atoms = [atom for res in mov_res for atom in res.get_atoms()]
    else:
        ref_atoms = list(ref_chain.get_atoms())
        mov_atoms = list(mov_chain.get_atoms())

    # Check if fixed and moving atoms lists have different size
    ref_atoms, mov_atoms = same_number_of_atoms(ref_atoms, mov_atoms)

    # Superposition between model and chain
    superimposer = Superimposer()
    superimposer.set_atoms(ref_atoms, mov_atoms)
    rmsd = superimposer.rms
    if ref_chain.xtra["type"] == "nuc":
        rmsd = rmsd / 2

    return rmsd, superimposer.rotran


def same_number_of_atoms(atom_chain1, atom_chain2):
    if len(atom_chain1) != len(atom_chain2):
        # Get the same number of atoms
        atoms_length = min(len(atom_chain1), len(atom_chain2))
        atom_chain2_idx, remainder = divmod((len(atom_chain2) - atoms_length), 2)
        atom_chain2 = atom_chain2[atom_chain2_idx:len(atom_chain2) - atom_chain2_idx - remainder]
        atom_chain1_idx, remainder = divmod((len(atom_chain1) - atoms_length), 2)
        atom_chain1 = atom_chain1[atom_chain1_idx:len(atom_chain1) - atom_chain1_idx - remainder]
    return atom_chain1, atom_chain2




def save_pdb(structure, out_dir, complex_name):
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(os.path.join(out_dir, "structures", complex_name + ".cif"))
