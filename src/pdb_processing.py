from Bio.PDB import PDBParser, NeighborSearch, Selection, Superimposer
from Bio.pairwise2 import align, format_alignment
from string import ascii_uppercase

processed_chains = list()  # List to store the ModelChains that are created
chain_to_model_chain = dict()  # Dictionary to map every chain in the input files to a ModelChain

letter_list = list(ascii_uppercase)  # List of A-Z to use as ModelChain custom ids
chain_ids = letter_list + [a + b for a in letter_list for b in letter_list]  # Extended list A-ZZ


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
        if is_interaction(structure, identity_threshold, ns_threshold, rmsd_threshold):
            initialize_model_chains(structure, identity_threshold, ns_threshold, rmsd_threshold)
            add_interactions(structure)
        else:
            # Warning
            pass

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



def is_interaction(structure, identity_threshold, ns_threshold, rmsd_threshold):
    chain_list = list(structure.get_chains())
    # Get atoms from each pair of chains
    for chain in chain_list:
        chain1 = chain_list[0]
        chain2 = chain_list[1]
        atoms1 = list(chain1.get_atoms())
        atoms2 = list(chain2.get_atoms())
        if len(chain_list) == 3:
            chain3 = chain_list[2]
            atoms3 = list(chain3.get_atoms())
            neighbors_search3 = NeighborSearch(atoms3)
            for atom in atoms1:
                close_atoms13 = neighbors_search3.search(atom.coord, float(ns_threshold))
                if len(close_atoms13) > 0:
                    return True
                else:
                    continue
                return False
            for atom in atoms2:
                close_atoms23 = neighbors_search3.search(atom.coord, float(ns_threshold))
                if len(close_atoms23) > 0:
                    return True
                else:
                    continue
                return False
        # Neighbor Search to find atoms within a given threshold
        neighbors_search1 = NeighborSearch(atoms1)
        for atom in atoms2:
            # List of atoms closer than ns_threshold
            close_atoms12 = neighbors_search1.search(atom.coord, float(ns_threshold))
            # Check that at least one pair of atoms are interacting
            if len(close_atoms12) > 0:
                return True
            else:
                continue
            return False



def initialize_model_chains(structure, identity_threshold, ns_threshold, rmsd_threshold):
    # Make a list of Biopython chain objects from the Structure
    chain_list = list(structure.get_chains())

    # For each PDB, go over its chains and create the ModelChains
    for chain in chain_list:
        # If this isn't the first PDB to be processed
        if len(processed_chains) > 0:
            # Check if there's any ModelChain that matches or not
            similar_chain_model = get_similar_chain_model(chain, identity_threshold, ns_threshold,
                                                          rmsd_threshold)
            if similar_chain_model is None:  # If there isn't, create a new ModelChain
                sequence = chain.xtra
                new_model_chain = ModelChain(chain, sequence)
                processed_chains.append(new_model_chain)
                chain_id = ":".join(chain.get_full_id()[0::2])
                chain_to_model_chain[chain_id] = new_model_chain
            else:  # If there is, map the chain to the existing ModelChain
                chain_id = ":".join(chain.get_full_id()[0::2])
                chain_to_model_chain[chain_id] = similar_chain_model
        # If it is the first PDB to be processed
        else:
            # Create the first ModelChain
            sequence = chain.xtra
            new_model_chain = ModelChain(chain, sequence)
            processed_chains.append(new_model_chain)
            chain_id = ":".join(chain.get_full_id()[0::2])
            chain_to_model_chain[chain_id] = new_model_chain



def get_similar_chain_model(chain, identity_threshold, ns_threshold, rmsd_threshold):
    # Get chain sequence from the xtra attribute
    chain_seq = chain.xtra
    # Compare it with each ModelChain that exists and return it if there's one that matches
    for model_chain in processed_chains:
        # Compare sequence similarity
        if align_with_modelchain(chain_seq, model_chain, identity_threshold):
            # If sequences are similar enough, superimpose the structures
            if superimpose_with_modelchain(chain, model_chain, rmsd_threshold):
                # If both are similar, return the ModelChain to which the chain will belong
                return model_chain
    # If there isn't any, return None
    return None

def align_with_modelchain(chain_seq, model_chain, identity_threshold):
    # Compare the chain_seq with the sequence of a ModelChain
    alignment = align.globalxx(chain_seq, model_chain.sequence, one_alignment_only=True)[0]
    longest_length = max(len(chain_seq), len(model_chain.sequence))
    if (alignment.score / longest_length) > float(identity_threshold):
        return True
    return False

def superimpose_with_modelchain(chain, model_chain, rmsd_threshold):
    # Calculate RMSD between the homologous chains
    chain_atoms = list(chain.get_atoms())
    model_atoms = list(model_chain.chain.get_atoms())
    # Check if fixed and moving atoms lists have different size
    if len(model_atoms) != len(chain_atoms):
        # Get the same amount of atoms
        atoms_length = min(len(model_atoms), len(chain_atoms))
        chain_atoms_idx = (len(chain_atoms) - atoms_length)/2
        chain_atoms = chain_atoms[chain_atoms_idx:len(chain_atoms)-chain_atoms_idx]
        model_atoms_idx = (len(model_atoms) - atoms_length)/2
        model_atoms = model_atoms[model_atoms_idx:len(model_atoms)-model_atoms_idx]
    # Superposition between model and chain
    super_imposer = Superimposer()
    super_imposer.set_atoms(model_atoms, chain_atoms)
    RMSD = super_imposer.rms
    if RMSD < rmsd_threshold:
        return True
    return False



def add_interactions(structure):
    chain_list = list(structure.get_chains())

    # As a first approximation, suppose that all chains in the PDB interact; then for each chain
    for chain in chain_list:
        chain_id = ":".join(chain.get_full_id()[0::2])
        # Get the corresponding ModelChain
        model_chain = chain_to_model_chain[chain_id]
        # Go over the rest of the chains in the file and add them to the ModelChain's interactions
        for interacting_chain in chain_list:
            if interacting_chain is not chain:
                interacting_chain_id = ":".join(interacting_chain.get_full_id()[0::2])
                interacting_model_chain = chain_to_model_chain[interacting_chain_id]
                model_chain.add_interaction(chain, interacting_model_chain, interacting_chain)
