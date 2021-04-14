import logging as logging
import sys
# from PDB_processing import processed_chains, chain_to_model_chain, model_chain_to_chains
from PDB_processing import chain_to_model_chain, model_chain_to_chains
from PDB_tools import get_chain_full_id



def parse_stoichiometry(stoichiometry_path, all_chains_by_pdb, all_prefixes, all_chains, stoich_dict):
    bad_format_msg = f"Error in stoichiometry file. The format should be one of the following:\n" \
                     "<Chain ID>:<number>\n" \
                     "<PDB ID>.<Chain ID>:<number>\n" \
                     "<Uniprot ID>:<number>"
    with open(stoichiometry_path, 'r') as stoich_file:
        lines = stoich_file.readlines()
        for line in lines:
            stoich_id = line.split(":")[0]
            num = line.split(":")[1]
            if "." in stoich_id:
                struct_id = stoich_id.split(".")[0]
                chain_id = stoich_id.split(".")[1]
                if struct_id is None or chain_id is None or num is None or \
                        not struct_id.isalnum() or not chain_id.isalpha() or \
                        not num.strip().isnumeric():
                    log.error(bad_format_msg)
                    sys.exit(1)

                if struct_id not in all_chains_by_pdb:
                    log.error(f"Error in stoichiometry file: structure {struct_id} is not present")
                    sys.exit(1)
                if chain_id not in all_chains_by_pdb[struct_id]:
                    log.error(f"Error in stoichiometry file: chain {chain_id} is not present"
                              f" in structure {struct_id}")
                    sys.exit(1)
            else:
                if stoich_id is None or num is None or \
                        not stoich_id.isalnum() or not num.strip().isnumeric():
                    log.error(bad_format_msg)
                    sys.exit(1)
                if stoich_id not in all_prefixes and stoich_id not in all_chains:
                    log.error(f"Error in stoichiometry file: chain/prefix {stoich_id} is not present in any"
                              " structure")
                    sys.exit(1)

            stoich_dict[stoich_id] = int(num)




# Check if the chain is present in the stoichiometry
def is_in_stoichiometry(stoich_dict, chain):
    if stoichiometry_ids_are_prefixes(stoich_dict):
        stoich_id = get_chain_full_id(chain).split(".")[0]
        return stoich_id in stoich_dict

    model_chain = chain_to_model_chain[get_chain_full_id(chain)]
    chains = model_chain_to_chains[model_chain.id]
    if stoichiometry_ids_are_struct_chain(stoich_dict):
        struct_id = get_chain_full_id(chain).split("_")[0]
        return any([struct_id + "." + c in stoich_dict for c in chains])

    return any([chain in stoich_dict for chain in chains])



# Get the chain id used in the stoichiometry to use as common chain id
def get_common_chain_id(stoich_dict, chain):
    if stoichiometry_ids_are_prefixes(stoich_dict):
        stoich_id = get_chain_full_id(chain).split(".")[0]
        if stoich_id in stoich_dict:
            return stoich_id
    else:
        if chain.get_id() in stoich_dict:
            return chain.get_id()

        model_chain = chain_to_model_chain[get_chain_full_id(chain)]
        chains = model_chain_to_chains[model_chain.id]
        if stoichiometry_ids_are_struct_chain(stoich_dict):
            struct_id = get_chain_full_id(chain).split("_")[0]
            return next(struct_id + "." + c for c in chains
                        if struct_id + "." + c in stoich_dict)

        return next(c for c in chains if c in stoich_dict)



def stoichiometry_ids_are_prefixes(stoich_dict):
    prefix_len = 6
    return any(
        len(stoich_id) == prefix_len and stoich_id.isalnum() for stoich_id in stoich_dict.keys())



def stoichiometry_ids_are_struct_chain(stoich_dict):
    return any("." in stoich_id for stoich_id in stoich_dict.keys())



def update_stoichiometry(stoich_dict, chain, current_stoich_dict, current_stoich_dict_by_prefix):
    if stoichiometry_ids_are_prefixes(stoich_dict):
        stoich_id = get_chain_full_id(chain).split(".")[0]
        current_stoich_dict_by_prefix[stoich_id] = current_stoich_dict_by_prefix.setdefault(
            stoich_id, dict())
        current_stoich_dict_by_prefix[stoich_id][
            chain.get_id()] = current_stoich_dict_by_prefix[stoich_id].setdefault(chain.get_id(),
                                                                                  0) + 1

        structure = chain.parent.parent
        structure_chain_ids = [c.get_id() for c in structure.get_chains()]
        if list(current_stoich_dict_by_prefix[stoich_id].keys()) == structure_chain_ids and all(
                [current_stoich_dict_by_prefix[stoich_id][c.get_id()] == 1 for c in
                 structure.get_chains()]):
            current_stoich_dict[stoich_id] = current_stoich_dict.setdefault(
                stoich_id, 0) + 1
            for chain_id in current_stoich_dict_by_prefix[stoich_id].keys():
                current_stoich_dict_by_prefix[stoich_id][chain_id] = 0
    else:
        common_chain_id = get_common_chain_id(stoich_dict, chain)
        current_stoich_dict[common_chain_id] = current_stoich_dict.setdefault(
            common_chain_id, 0) + 1
