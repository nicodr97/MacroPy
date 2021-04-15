import logging as log
from Bio.PDB import NeighborSearch
from string import ascii_uppercase
from PDB_tools import align, get_chain_full_id, superimpose

processed_chains = list()  # List to store the ModelChains that are created
chain_to_model_chain = dict()  # Dictionary to map every chain in the input files to a ModelChain
model_chain_to_chains = dict()  # Dictionary to map every ModelChain to its homolog/similar chains

letter_list = list(ascii_uppercase)  # List of A-Z to use as ModelChain custom ids
modelchain_ids = letter_list + [a + b for a in letter_list for b in letter_list]  # Large list A-ZZ


# Definition of the class ModelChain
class ModelChain:

    def __init__(self, chain):
        self.id = modelchain_ids.pop(0)
        self.chain = chain
        self.sequence = chain.xtra["seq"]
        self.interactions = list()

    def add_interaction(self, homolog_chain, interacting_chain):
        self.interactions.append((homolog_chain, interacting_chain))


def process_pdbs(pdb_chains, identity_threshold, ns_threshold, rmsd_threshold):
    """Process PDB's and store their information in ModelChain objects"""
    log.info("Processing chains")
    # Go over all the PDB chains in pdb_chains, processing from longer to shorter
    for chain in sorted(pdb_chains, key=lambda x: (len(x.xtra["seq"]), x.xtra["full_id"]),
                        reverse=True):
        initialize_model_chains(chain, identity_threshold, rmsd_threshold)
        add_interactions(chain, ns_threshold)


def initialize_model_chains(chain, identity_threshold, rmsd_threshold):
    """Create ModelChain objects with the processed chains from the PDB files"""

    # If this isn't the first chain to be processed
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
    """
    Get ModelChain similar to the given chain based on the sequence identity
    and the RMSD between them
    """

    for model_chain in filter(lambda x: x.chain.xtra["type"] == chain.xtra["type"],
                              processed_chains):
        # Align the sequences to compare the alignment score %
        # with the identity threshold
        score, positions = next(align(chain, model_chain.sequence))
        if score > identity_threshold:
            # Superimpose them to compare the RMSD with its threshold
            rmsd, rotran = superimpose(chain, model_chain.chain, positions)
            if rmsd < rmsd_threshold:
                return model_chain
    # If there isn't any, return None
    return None


def add_interactions(chain, ns_threshold):
    """Add interactions to ModelChain"""

    # Get all the chains of the PDB from which the processed chain comes
    chain_list = [int_chain for int_chain in chain.parent.parent.get_chains()
                  if int_chain is not chain]
    chain_id = get_chain_full_id(chain)

    # Get the corresponding ModelChain of the chain
    model_chain = chain_to_model_chain[chain_id]
    # Go over the rest of the chains in the PDB and add them to the ModelChain's interactions
    for interacting_chain in chain_list:
        if are_interacting(chain, interacting_chain, ns_threshold):
            model_chain.add_interaction(chain, interacting_chain)


def are_interacting(chain1, chain2, ns_threshold):
    """
    Check if two chains are interacting,
    i.e. if at least one pair of atoms are closer than a threshold
    """

    # Get atoms from each chain
    atoms1 = list(chain1.get_atoms())
    atoms2 = list(chain2.get_atoms())

    neighbors_search1 = NeighborSearch(atoms1)
    for atom in atoms2:
        # List of atoms closer than ns_threshold
        close_atoms1_2 = neighbors_search1.search(atom.coord, ns_threshold)
        # Check that at least one pair of atoms are interacting
        if len(close_atoms1_2) > 0:
            return True
    return False
