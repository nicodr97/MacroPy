
from Bio.PDB import NeighborSearch, Superimposer
# from Bio.pairwise2 import align
from Bio.Align import PairwiseAligner
from string import ascii_uppercase

get_chain_full_id = lambda chain: ":".join(chain.get_full_id()[0::2])
processed_chains = list()  # List to store the ModelChains that are created
chain_to_model_chain = dict()  # Dictionary to map every chain in the input files to a ModelChain

letter_list = list(ascii_uppercase)  # List of A-Z to use as ModelChain custom ids
modelchain_ids = letter_list + [a + b for a in letter_list for b in letter_list]  # Large list A-ZZ


# Definition of the class ModelChain
class ModelChain:

    def __init__(self, chain, sequence):
        self.id = modelchain_ids.pop(0)
        self.chain = chain
        self.sequence = sequence
        self.interactions = list()

    def add_interaction(self, homolog_chain, interacting_model_chain, interacting_chain):
        self.interactions.append((homolog_chain, interacting_model_chain, interacting_chain))


# Process PDBs and store their information in ModelChain objects
def process_pdbs(pdb_dict, identity_threshold, ns_threshold, rmsd_threshold):
    # Go over all the PDBs in the pdb_dict
    for pdb_id, structure in pdb_dict.items():
        initialize_model_chains(structure, identity_threshold, rmsd_threshold)
        add_interactions(structure, ns_threshold)

    ######## Loop for checking the result
    print("Processed chains: \n\n")
    print(chain_to_model_chain)
    for model_chain in processed_chains:
        print(model_chain.id)
        print(model_chain.sequence)
        for chain_interactions in model_chain.interactions:
            chain1 = chain_interactions[0]
            interacting_model_chain = chain_interactions[1]
            chain2 = chain_interactions[2]
            chain1_structure = chain1.parent.parent
            chain2_structure = chain2.parent.parent

            print(f"Interaction of {chain1.get_id()} of structure {chain1_structure.get_id()} and "
                  f"{chain2.get_id()} of structure {chain2_structure.get_id()} and model chain "
                  f"{interacting_model_chain.id}")


def initialize_model_chains(structure, identity_threshold, rmsd_threshold):

    # Make a list of Biopython chain objects from the Structure
    chain_list = list(structure.get_chains())

    # For each PDB, go over its chains and create the ModelChains
    for chain in chain_list:
        # if chain.xtra["type"] == "prot":
        # If this isn't the first PDB to be processed
        if len(processed_chains) > 0:
            # Check if there's any ModelChain that matches or not
            similar_chain_model = get_similar_chain_model(chain, identity_threshold, rmsd_threshold)
            if similar_chain_model is None:  # If there isn't, create a new ModelChain
                sequence = chain.xtra["seq"]
                new_model_chain = ModelChain(chain, sequence)
                processed_chains.append(new_model_chain)
                chain_id = get_chain_full_id(chain)
                chain_to_model_chain[chain_id] = new_model_chain
            else:  # If there is, map the chain to the existing ModelChain
                if len(chain.xtra["seq"]) > len(similar_chain_model.sequence):
                    similar_chain_model.sequence = chain.xtra["seq"]
                    similar_chain_model.chain = chain
                chain_id = get_chain_full_id(chain)
                chain_to_model_chain[chain_id] = similar_chain_model
        # If it is the first PDB to be processed
        else:
            # Create the first ModelChain
            sequence = chain.xtra["seq"]
            new_model_chain = ModelChain(chain, sequence)
            processed_chains.append(new_model_chain)
            chain_id = get_chain_full_id(chain)
            chain_to_model_chain[chain_id] = new_model_chain

        # elif chain.xtra["type"] == "nuc":
        #     chain_seq = chain.xtra["seq"]
        #     print(f"checking {get_chain_full_id(chain)}: {chain_seq}")
        #     if len(processed_chains) > 0:
        #         for modelchain_obj in (modchain for modchain in processed_chains if modchain.chain.xtra["type"] == "nuc"):
        #             print(f"against modelchain: {modelchain_obj.sequence}")
        #             if chain_seq[1:-1] in modelchain_obj.sequence:
        #                 print(f"chain_seq {chain_seq} is inside modelchain.sequence {modelchain_obj.sequence}\n\n")
        #                 chain_id = get_chain_full_id(chain)
        #                 chain_to_model_chain[chain_id] = modelchain_obj
        #                 break
        #             elif modelchain_obj.sequence[1:-1] in chain_seq:
        #                 print(f"modelchain.sequence {modelchain_obj.sequence} is inside chain_seq {chain_seq}\n\n")
        #                 modelchain_obj.chain = chain
        #                 modelchain_obj.sequence = chain_seq
        #                 chain_id = get_chain_full_id(chain)
        #                 chain_to_model_chain[chain_id] = modelchain_obj
        #                 break
        #     print("They are different\n\n")
        #     new_model_chain = ModelChain(chain, chain_seq)
        #     processed_chains.append(new_model_chain)
        #     chain_id = get_chain_full_id(chain)
        #     chain_to_model_chain[chain_id] = new_model_chain



def get_similar_chain_model(chain, identity_threshold, rmsd_threshold):
    print(f"\nget similar model chain of: {chain.xtra['seq']}")
    # Compare it with each ModelChain that exists and return it if there's one that matches
    for model_chain in processed_chains:
        # If sequences are similar enough, superimpose the structures
        are_aligned, positions = align_with_modelchain(chain, model_chain, identity_threshold)
        print(f"comparing with {model_chain.sequence}: is aligned {are_aligned} ({positions})")
        if are_aligned and superimpose_with_modelchain(chain, model_chain, positions, rmsd_threshold):
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

    # alignment = align.localms(chain_seq, model_chain.sequence, 1, -1, -1, -1, one_alignment_only=True)[0]
    # longest_length = max(len(chain_seq), len(model_chain.sequence))
    # return (alignment.score / longest_length) > float(identity_threshold)


def superimpose_with_modelchain(chain, model_chain, positions, rmsd_threshold):
    chain_idx = positions[0][0]
    chain_res = list(chain.get_residues())[chain_idx[0]:chain_idx[1]]
    chain_atoms = [atom for res in chain_res for atom in list(res.get_atoms())]
    model_idx = positions[1][0]
    model_res = list(model_chain.chain.get_residues())[model_idx[0]:model_idx[1]]
    model_atoms = [atom for res in model_res for atom in list(res.get_atoms())]
    # Calculate RMSD between the homologous chains
    # chain_atoms = list(chain.get_atoms())
    # model_atoms = list(model_chain.chain.get_atoms())
    # Check if fixed and moving atoms lists have different size
    if len(model_atoms) != len(chain_atoms):
        # Get the same amount of atoms
        atoms_length = min(len(model_atoms), len(chain_atoms))
        chain_atoms_idx, remainder = divmod((len(chain_atoms) - atoms_length), 2)
        chain_atoms = chain_atoms[chain_atoms_idx:len(chain_atoms) - chain_atoms_idx - remainder]
        model_atoms_idx, remainder = divmod((len(model_atoms) - atoms_length), 2)
        model_atoms = model_atoms[model_atoms_idx:len(model_atoms) - model_atoms_idx - remainder]

    # Superposition between model and chain
    super_imposer = Superimposer()
    super_imposer.set_atoms(model_atoms, chain_atoms)
    rmsd = super_imposer.rms
    if chain.xtra["type"] == "nuc":
        rmsd = rmsd/2
    print(f"{rmsd}")

    return rmsd < float(rmsd_threshold)


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
                interacting_model_chain.add_interaction(interacting_chain, model_chain, chain)

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
