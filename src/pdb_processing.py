from Bio.PDB import PDBParser, NeighborSearch, Selection, Polypeptide
from Bio.pairwise2 import align, format_alignment

# Dictionary to store the ModelChains that are created
processed_chains = dict()

# Definition of the class ModelChain
class ModelChain(object):

    def __init__(self):
        self.chains = list()
        # self.interactions = dict()

    def add_chain(self, chain_obj):
        self.chains.append(chain_obj)
    # def add_interaction(self, modelchain_obj, chain_obj):
    #     self.interactions[modelchain_obj] = chain_obj

    def get_sequences(self):
        for chain in self.chains:
            yield Polypeptide.PPBuilder().build_peptides(chain)[0].get_sequence()


# Process PDBs and store their information in ModelChain objects
def process_pdbs(pdb_dict, id_thres, rmsd_thres):
    i=1

    # Go over all the PDBs in the pdb_dict
    for pdb_id, structure_obj in pdb_dict.items():

        ####### Check that they are interacting: NegihborSearch. If they are, proceed with:

        # For each PDB, go over its chains
        chain_list = structure_obj.get_chains()
        for chain in structure_obj.get_chains():
            # If this isn't the first PDB to be processed
            if len(processed_chains) > 0:
                # Compare each chain with each of the already processed ModelChains
                for modelchain_obj in processed_chains.values():
                    # Get sequence of the chain with Polypeptide.PPBuilder
                    chain_seq = Polypeptide.PPBuilder().build_peptides(chain)[0].get_sequence()
                    if compare_with_modelchain(chain_seq, modelchain_obj, id_thres, rmsd_thres):
                        # If it meets the requirements, add it to the ModelChain
                        modelchain_obj.add_chain(chain)
                        break
                    else:
                        # If it doesn't meet the requirements, create a new ModelChain and add it
                        processed_chains[i] = ModelChain()
                        processed_chains[i].add_chain(chain)
                        i+=1
                        break

            # If it is the first PDB to be processed
            else:
                # Create the first ModelChain
                processed_chains[i] = ModelChain()
                processed_chains[i].add_chain(chain)
                i+=1

    ######### Loop for checking the result
    for modchain in processed_chains.values():
        print(modchain)
        for chain in modchain.chains:
            print(chain.get_id())


def compare_with_modelchain(chain_seq, modelchain_obj, id_thres, rmsd_thres):
    # To compare the chain_seq with every sequence in a ModelChain object
    for seq in modelchain_obj.get_sequences():
        alignment = align.globalxx(chain_seq, seq, one_alignment_only=True)[0]
        if (alignment.score / len(chain_seq)) > id_thres:  # or max(len(seq1), len(seq2))?

        ####### RMSD con args.RMSD-threshold

            continue
        else:
            return False
    return True
