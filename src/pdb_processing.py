from Bio.PDB import PDBParser, NeighborSearch, Selection, Polypeptide
from Bio.pairwise2 import align, format_alignment
from string import ascii_uppercase

processed_chains = list()  # Dictionary to store the ModelChains that are created
chain_ids = list(ascii_uppercase)
chain_to_model_chain = dict()


# Definition of the class ModelChain
class ModelChain:

    def __init__(self, chain):
        self.id = chain_ids.pop(0)
        self.chain = chain
        self.interactions = list()

    def add_interaction(self, homolog_chain, interacting_model_chain, interacting_chain):
        self.interactions.append((homolog_chain, interacting_model_chain, interacting_chain))

    def get_sequence(self):
        return Polypeptide.PPBuilder().build_peptides(self.chain)[0].get_sequence()


# Process PDBs and store their information in ModelChain objects
def process_pdbs(pdb_dict, identity_threshold, rmsd_threshold):
    # Go over all the PDBs in the pdb_dict
    for pdb_id, structure in pdb_dict.items():
        ####### Check that they are interacting: NegihborSearch. If they are, proceed with:
        chain_list = list(structure.get_chains())

        initialize_model_chains(chain_list, identity_threshold, rmsd_threshold)
        add_interactions(chain_list)

    print("Processed chains: \n\n")
    ######### Loop for checking the result
    for model_chain in processed_chains:
        print(model_chain.id)
        for chain_interactions in model_chain.interactions:
            chain1 = chain_interactions[0]
            interacting_model_chain = chain_interactions[1]
            chain2 = chain_interactions[2]
            chain1_structure = Selection.unfold_entities(chain1, "S")[0]
            chain2_structure = Selection.unfold_entities(chain2, "S")[0]

            print(f"Interaction of {chain1.get_id()} of structure {chain1_structure.get_id()} and "
                  f"{chain2.get_id()} of structure {chain2_structure.get_id()} and model chain "
                  f"{interacting_model_chain.id}")


def compare_with_model_chain(chain_seq, model_chain, identity_threshold, rmsd_threshold):
    # Compare the chain_seq with the sequence of a ModelChain
    alignment = align.globalxx(chain_seq, model_chain.get_sequence(), one_alignment_only=True)[0]
    if (alignment.score / len(chain_seq)) > identity_threshold:  # or max(len(seq1), len(seq2))?
        return True

        ####### RMSD con args.RMSD-threshold

    return False


def chain_model_exists(chain, identity_threshold, rmsd_threshold):
    chain_seq = Polypeptide.PPBuilder().build_peptides(chain)[0].get_sequence()
    # Compare each chain with each of the already processed ModelChains
    return any([compare_with_model_chain(chain_seq, model_chain, identity_threshold,
                                         rmsd_threshold) for model_chain in processed_chains])


def get_similar_chain_model(chain, identity_threshold, rmsd_threshold):
    chain_seq = Polypeptide.PPBuilder().build_peptides(chain)[0].get_sequence()
    for model_chain in processed_chains:
        if compare_with_model_chain(chain_seq, model_chain, identity_threshold,
                                    rmsd_threshold):
            return model_chain

    return None


def initialize_model_chains(chain_list, identity_threshold, rmsd_threshold):
    # For each PDB, go over its chains and create the ModelChains
    for chain in chain_list:
        # If this isn't the first PDB to be processed
        if len(processed_chains) > 0:
            similar_chain_model = get_similar_chain_model(chain, identity_threshold, rmsd_threshold)
            if similar_chain_model is None:
                new_model_chain = ModelChain(chain)
                processed_chains.append(new_model_chain)
                chain_id = get_chain_global_id(chain)
                chain_to_model_chain[chain_id] = new_model_chain
            else:
                chain_id = get_chain_global_id(chain)
                chain_to_model_chain[chain_id] = similar_chain_model
        # If it is the first PDB to be processed
        else:
            # Create the first ModelChain
            new_model_chain = ModelChain(chain)
            processed_chains.append(new_model_chain)
            chain_id = get_chain_global_id(chain)
            chain_to_model_chain[chain_id] = new_model_chain


def get_chain_global_id(chain):
    chain_structure = Selection.unfold_entities(chain, "S")[0]
    return chain_structure.get_id() + "_" + chain.get_id()


def add_interactions(chain_list):
    # As a first approximation, suppose that all chains in the PDB interact
    for chain in chain_list:
        chain_id = get_chain_global_id(chain)
        model_chain = chain_to_model_chain[chain_id]
        for interacting_chain in chain_list:
            if interacting_chain is not chain:
                interacting_chain_id = get_chain_global_id(interacting_chain)
                interacting_model_chain = chain_to_model_chain[interacting_chain_id]
                model_chain.add_interaction(chain, interacting_model_chain, interacting_chain)
