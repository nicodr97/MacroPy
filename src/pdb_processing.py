from Bio.PDB import PDBParser, NeighborSearch, Selection, Polypeptide
from Bio.pairwise2 import align, format_alignment
from string import ascii_uppercase

processed_chains = list()  # List to store the ModelChains that are created
chain_to_model_chain = dict() # Dictionary to map every chain in the input files to a ModelChain

letter_list = list(ascii_uppercase) # List of A-Z to use as ModelChain custom ids
chain_ids = letter_list + [a+b for a in letter_list for b in letter_list] # Extended list A-ZZ


# Definition of the class ModelChain
class ModelChain:

    def __init__(self, chain, sequence):
        self.id = chain_ids.pop(0)
        self.chain = chain
        self.sequence = sequence
        self.interactions = list()

    def add_interaction(self, homolog_chain, interacting_model_chain, interacting_chain):
        self.interactions.append((homolog_chain, interacting_model_chain, interacting_chain))



# Process PDBs and store their information in ModelChain objects
def process_pdbs(pdb_dict, identity_threshold, ns_threshold, rmsd_threshold):
    # Go over all the PDBs in the pdb_dict
    for pdb_id, structure in pdb_dict.items():

        # Get atoms from each pair of chains
        chain1, chain2 = list(structure.get_chains())
        atoms1 = list(chain1.get_atoms())
        atoms2 = list(chain2.get_atoms())

        # Neighbor Search to find atoms within a given threshold
        neighbors_search = NeighborSearch(atoms1)
        for atom in atoms2:
            # List of atoms closer than ns_threshold
            close_atoms = neighbors_search.search(atom.coord, float(ns_threshold))

        # Check that at least one pair of atoms are interacting
        if len(close_atoms) > 0:
            initialize_model_chains(structure, identity_threshold, ns_threshold, rmsd_threshold)
            add_interactions(structure)
        
    ######### Loop for checking the result
    print("Processed chains: \n\n")
    print(chain_to_model_chain)
    for model_chain in processed_chains:
        print(model_chain.id)
        print(model_chain.sequence)
        for chain_interactions in model_chain.interactions:
            chain1 = chain_interactions[0]
            interacting_model_chain = chain_interactions[1]
            chain2 = chain_interactions[2]
            chain1_structure = Selection.unfold_entities(chain1, "S")[0]
            chain2_structure = Selection.unfold_entities(chain2, "S")[0]

            print(f"Interaction of {chain1.get_id()} of structure {chain1_structure.get_id()} and "
                  f"{chain2.get_id()} of structure {chain2_structure.get_id()} and model chain "
                  f"{interacting_model_chain.id}")



# def chain_model_exists(chain, identity_threshold, rmsd_threshold):
#     chain_seq = Polypeptide.PPBuilder().build_peptides(chain)[0].sequence
#     # Compare each chain with each of the already processed ModelChains
#     return any([compare_with_model_chain(chain_seq, model_chain, identity_threshold,
#                                          rmsd_threshold) for model_chain in processed_chains])



def initialize_model_chains(structure, identity_threshold, ns_threshold, rmsd_threshold):
    # Make a list of Biopython chain objects from the Structure
    chain_list = list(structure.get_chains())

    # For each PDB, go over its chains and create the ModelChains
    for chain in chain_list:
        chain.xtra = structure.xtra # Transfer the biskit's sequences from Structure to Chain
        # If this isn't the first PDB to be processed
        if len(processed_chains) > 0:
            # Check if there's any ModelChain that matches or not
            similar_chain_model = get_similar_chain_model(chain, identity_threshold, ns_threshold, rmsd_threshold)
            if similar_chain_model is None: # If there isn't, create a new ModelChain
                sequence = chain.xtra[chain.get_id()]
                new_model_chain = ModelChain(chain, sequence)
                processed_chains.append(new_model_chain)
                chain_id = ":".join( chain.get_full_id()[0::2] )
                chain_to_model_chain[chain_id] = new_model_chain
            else: # If there is, map the chain to the existing ModelChain
                chain_id = ":".join( chain.get_full_id()[0::2] )
                chain_to_model_chain[chain_id] = similar_chain_model
        # If it is the first PDB to be processed
        else:
            # Create the first ModelChain
            sequence = chain.xtra[chain.get_id()]
            new_model_chain = ModelChain(chain, sequence)
            processed_chains.append(new_model_chain)
            chain_id = ":".join( chain.get_full_id()[0::2] )
            chain_to_model_chain[chain_id] = new_model_chain



def get_similar_chain_model(chain, identity_threshold, ns_threshold, rmsd_threshold):
    # Get biskit's chain sequence from the xtra attribute
    chain_seq = chain.xtra[chain.get_id()]

    # Compare it with each ModelChain that exists and return it if there's one that matches
    for model_chain in processed_chains:
        if compare_with_model_chain(chain_seq, model_chain, identity_threshold, 
                                    ns_threshold, rmsd_threshold):
            return model_chain

    # If there isn't any, return None
    return None


def compare_with_model_chain(chain_seq, model_chain, identity_threshold, ns_threshold, rmsd_threshold):
    # Compare the chain_seq with the sequence of a ModelChain
    alignment = align.globalxx(chain_seq, model_chain.sequence, one_alignment_only=True)[0]
    if (alignment.score / len(chain_seq)) > float(identity_threshold):  # or max(len(seq1), len(seq2))?
        return True

        ####### RMSD con args.RMSD-threshold

    return False



# def get_chain_global_id(chain):
#     chain_structure = Selection.unfold_entities(chain, "S")[0]
#     return chain_structure.get_id() + "_" + chain.get_id()





def add_interactions(structure):
    chain_list = list(structure.get_chains())

    # As a first approximation, suppose that all chains in the PDB interact; then for each chain
    for chain in chain_list:
        chain.xtra = structure.xtra # Transfer the biskit's sequences from Structure to Chain
        chain_id = ":".join( chain.get_full_id()[0::2] )
        # Get the corresponding ModelChain
        model_chain = chain_to_model_chain[chain_id]
        # Go over the rest of the chains in the file and add them to the ModelChain's interactions
        for interacting_chain in chain_list:
            if interacting_chain is not chain:
                interacting_chain_id = ":".join( interacting_chain.get_full_id()[0::2] )
                interacting_model_chain = chain_to_model_chain[interacting_chain_id]
                model_chain.add_interaction(chain, interacting_model_chain, interacting_chain)
