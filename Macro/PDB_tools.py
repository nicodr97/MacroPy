import os
import logging as log
from Bio.Align import PairwiseAligner
from Bio.PDB import Superimposer, MMCIFIO, PDBIO


def get_chain_full_id(chain):
    return ":".join(chain.get_full_id()[0::2])


def align(chain, seq):
    """Align the sequence of a chain object with a sequence string"""
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
        yield (i.score / length), i.aligned


def superimpose(ref_chain, mov_chain, alignment_positions=None):
    """
    Superimpose the chain to the reference chain,
    only by the aligned positions if they are provided
    """
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
    """Return two chains of atoms of equals size"""

    if len(atom_chain1) != len(atom_chain2):
        # Get the most centered atoms of each chain, so we don't take much more atoms of any
        # extreme of the chains
        atoms_length = min(len(atom_chain1), len(atom_chain2))
        atom_chain2_idx, remainder = divmod((len(atom_chain2) - atoms_length), 2)
        atom_chain2 = atom_chain2[atom_chain2_idx:len(atom_chain2) - atom_chain2_idx - remainder]
        atom_chain1_idx, remainder = divmod((len(atom_chain1) - atoms_length), 2)
        atom_chain1 = atom_chain1[atom_chain1_idx:len(atom_chain1) - atom_chain1_idx - remainder]
    return atom_chain1, atom_chain2


def save_structure(structure, out_dir, complex_name, save_pdb):
    """Save final complex"""

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(os.path.join(out_dir, "structures", complex_name + ".cif"))
    log.info(f"Complex saved as {complex_name}.cif")

    if save_pdb:
        io = PDBIO()
        chains_to_delete = [chain for chain in structure.get_chains()][26:]
        for chain in chains_to_delete:
            structure[0].detach_child(chain.get_id())
        io.set_structure(structure)
        io.save(os.path.join(out_dir, "structures", complex_name + ".pdb"))
        log.info(f"Complex saved as {complex_name}.pdb")


def minimize(out_dir, complex_name, max_steps):
    """Fill the missing atoms of the complex file and minimize it"""

    # Import only if the function is called
    from modeller import Environ, Selection
    from modeller import log as mlog
    from modeller.scripts import complete_pdb
    from modeller.optimizers import ConjugateGradients, actions

    # Set environment variables and import default libraries
    mlog.level(0)
    env = Environ(0)
    env.io.atom_files_directory = ['../atom_files']
    env.edat.dynamic_sphere = True
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Read the structure file, fill the missing atoms and create a selection (all)
    pdb_file = os.path.join(out_dir, "structures", complex_name + ".cif")
    log.info(f"Filling mising atoms of {complex_name}.cif prior to minimization")
    mdl = complete_pdb(env, pdb_file)
    atmsel = Selection(mdl)

    # Generate the minimization restraints
    mdl.restraints.make(atmsel, restraint_type='stereo', spline_on_site=False)

    # Create a ConjugateGradients object and perform the minimization
    cg = ConjugateGradients(output='REPORT')
    with open(os.path.join(out_dir, "analysis", "minimization" + ".log"), 'w') as trcfil:
        # Run until convergence or until max_steps and print progression to the output file
        # every 10 steps
        log.info(f"Running Conjugate Gradients minimization")
        cg.optimize(atmsel, min_atom_shift=0.01, max_iterations=max_steps,
                    actions=actions.Trace(10, trcfil))

    # Output the minimized structure
    mdl.write(file=os.path.join(out_dir, 'structures', complex_name + '_minimized.cif'),
              model_format="MMCIF")
    log.info(f"Minimized complex saved as {complex_name}_minimized.cif")
