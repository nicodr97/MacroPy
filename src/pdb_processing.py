import logging as log
from Bio.PDB import NeighborSearch, Superimposer
from Bio.Align import PairwiseAligner
from string import ascii_letters


processed_chains = list()  # List to store the ModelChains that are created
chain_to_model_chain = dict()  # Dictionary to map every chain in the input files to a ModelChain
model_chain_to_chains = dict()
modelchain_ids = list(ascii_letters)  # List of A-Z and a-z to use as ModelChain custom ids


# Definition of the class ModelChain
class ModelChain:

    def __init__(self, chain):
        self.id = modelchain_ids.pop(0)
        self.chain = chain
        self.sequence = chain.xtra["seq"]
        self.interactions = list()

    def add_interaction(self, homolog_chain, interacting_model_chain, interacting_chain):
        self.interactions.append((homolog_chain, interacting_model_chain, interacting_chain))


# Process PDBs and store their information in ModelChain objects
def process_pdbs(pdb_dict, identity_threshold, ns_threshold, rmsd_threshold):
    # Go over all the PDBs in the pdb_dict
    for pdb_id, structure in pdb_dict.items():
        initialize_model_chains(structure, identity_threshold, rmsd_threshold)
        add_interactions(structure, ns_threshold)

    # Loop for checking the result
    log.info("Processed chains: \n\n")
    log.info(chain_to_model_chain)
    for model_chain in processed_chains:
        log.info(model_chain.id)
        log.info(model_chain.sequence)
        for chain_interactions in model_chain.interactions:
            chain1 = chain_interactions[0]
            interacting_model_chain = chain_interactions[1]
            chain2 = chain_interactions[2]
            chain1_structure = chain1.parent.parent
            chain2_structure = chain2.parent.parent

            log.info(f"Interaction of {chain1.get_id()} of structure {chain1_structure.get_id()}"
                     f" and {chain2.get_id()} of structure {chain2_structure.get_id()} and model"
                     f" chain {interacting_model_chain.id}")


def initialize_model_chains(structure, identity_threshold, rmsd_threshold):
    # For each PDB, go over its chains and create the ModelChains
    for chain in structure.get_chains():
        # If this isn't the first PDB to be processed
        if len(processed_chains) > 0:
            # Check if there's any ModelChain that matches or not
            similar_chain_model = get_similar_chain_model(chain, identity_threshold, rmsd_threshold)
            if similar_chain_model is None:  # If there isn't, create a new ModelChain
                new_model_chain = ModelChain(chain)
                processed_chains.append(new_model_chain)
                chain_id = get_chain_full_id(chain)
                chain_to_model_chain[chain_id] = new_model_chain
                model_chain_to_chains[new_model_chain.id] = set(chain.get_id())
            else:  # If there is, map the chain to the existing ModelChain
                if len(chain.xtra["seq"]) > len(similar_chain_model.sequence):
                    similar_chain_model.sequence = chain.xtra["seq"]
                    similar_chain_model.chain = chain
                chain_id = get_chain_full_id(chain)
                chain_to_model_chain[chain_id] = similar_chain_model
                model_chain_to_chains[similar_chain_model.id].add(chain.get_id())
        # If it is the first PDB to be processed
        else:
            # Create the first ModelChain
            new_model_chain = ModelChain(chain)
            processed_chains.append(new_model_chain)
            chain_id = get_chain_full_id(chain)
            chain_to_model_chain[chain_id] = new_model_chain
            model_chain_to_chains[new_model_chain.id] = set(chain.get_id())


def get_similar_chain_model(chain, identity_threshold, rmsd_threshold):
    log.info(f"\nget similar model chain of: {chain.xtra['seq']}")
    # Compare it with each ModelChain that exists and return it if there's one that matches
    for model_chain in processed_chains:
        # If sequences are similar enough, superimpose the structures
        are_aligned, positions = align_with_modelchain(chain, model_chain, identity_threshold)
        log.info(f"comparing with {model_chain.sequence}: is aligned {are_aligned} ({positions})")
        if are_aligned and superimpose_with_modelchain(chain, model_chain, positions,
                                                       rmsd_threshold):
            # If both are similar, return the ModelChain to which the chain will belong
            return model_chain
    # If there isn't any, return None
    return None


def align_with_modelchain(chain, model_chain, identity_threshold):
    # Get chain sequence from the xtra attribute
    chain_seq = chain.xtra["seq"]
    # Compare the chain_seq with the sequence of a ModelChain
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    alignment = aligner.align(chain_seq, model_chain.sequence)[0]

    lengths = (len(chain_seq), len(model_chain.sequence))
    length = max(lengths) if chain.xtra["type"] == "prot" else min(lengths)
    are_aligned = (alignment.score / length) > float(identity_threshold)
    return are_aligned, alignment.aligned


def superimpose_with_modelchain(chain, model_chain, alignment_positions, rmsd_threshold):
    chain_idx = alignment_positions[0][0]
    start_res = chain_idx[0]
    end_res = chain_idx[1]
    chain_res = list(chain.get_residues())[start_res:end_res]
    chain_atoms = [atom for res in chain_res for atom in res.get_atoms()]

    model_idx = alignment_positions[1][0]
    start_res = model_idx[0]
    end_res = model_idx[1]
    model_res = list(model_chain.chain.get_residues())[start_res:end_res]
    model_atoms = [atom for res in model_res for atom in res.get_atoms()]

    # Calculate RMSD between the homologous chains
    # Check if fixed and moving atoms lists have different size
    model_atoms, chain_atoms = get_atom_chains_with_same_length(model_atoms, chain_atoms)

    # Superposition between model and chain
    super_imposer = Superimposer()
    super_imposer.set_atoms(model_atoms, chain_atoms)
    rmsd = super_imposer.rms
    if chain.xtra["type"] == "nuc":
        rmsd = rmsd / 2
    log.info(f"RMSD = {rmsd}")

    return rmsd < float(rmsd_threshold)


def get_atom_chains_with_same_length(atom_chain1, atom_chain2):
    if len(atom_chain1) != len(atom_chain2):
        # Get the same amount of atoms
        atoms_length = min(len(atom_chain1), len(atom_chain2))
        atom_chain2_idx, remainder = divmod((len(atom_chain2) - atoms_length), 2)
        atom_chain2 = atom_chain2[atom_chain2_idx:len(atom_chain2) - atom_chain2_idx - remainder]
        atom_chain1_idx, remainder = divmod((len(atom_chain1) - atoms_length), 2)
        atom_chain1 = atom_chain1[atom_chain1_idx:len(atom_chain1) - atom_chain1_idx - remainder]

    return atom_chain1, atom_chain2


def add_interactions(structure, ns_threshold):
    chain_list = list(structure.get_chains())

    for i in range(len(chain_list)):
        chain = chain_list[i]
        chain_id = get_chain_full_id(chain)
        # Get the corresponding ModelChain
        model_chain = chain_to_model_chain[chain_id]
        # Go over the rest of the chains in the file and add them to the ModelChain's interactions
        for interacting_chain in chain_list[i + 1:]:
            if are_interacting(chain, interacting_chain, ns_threshold):
                interacting_chain_id = get_chain_full_id(interacting_chain)
                interacting_model_chain = chain_to_model_chain[interacting_chain_id]
                model_chain.add_interaction(chain, interacting_model_chain, interacting_chain)
                if model_chain is not interacting_model_chain:
                    interacting_model_chain.add_interaction(interacting_chain, model_chain, chain)


def get_chain_full_id(chain):
    return ":".join(chain.get_full_id()[0::2])


def are_interacting(chain1, chain2, ns_threshold):
    # Get atoms from each pair of chains
    atoms1 = list(chain1.get_atoms())
    atoms2 = list(chain2.get_atoms())

    neighbors_search1 = NeighborSearch(atoms1)
    for atom in atoms2:
        # List of atoms closer than ns_threshold
        close_atoms1_2 = neighbors_search1.search(atom.coord, float(ns_threshold))
        # Check that at least one pair of atoms are interacting
        if len(close_atoms1_2) > 0:
            return True
    return False
