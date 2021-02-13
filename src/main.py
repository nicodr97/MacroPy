import argparse
import os
import sys
from Bio.PDB import PDBParser
import gzip


def parse_input_directory(path):
    if not os.path.isdir(path):
        print("Path is not a directory")
        return 1

    file_names = [f for f in os.listdir(path) if
                  os.path.isfile(os.path.join(path, f))]
    for file_name in file_names:
        file_name_parts = file_name.split(sep=".")
        error_msg = f"Error in file '{file_name}' . "
        wrong_naming_msg = "Wrong file naming structure. It should be as " \
                           "followed: <name>_<chain1>_<chain2>.pdb(.gz)"

        if ((len(file_name_parts) != 2 and len(file_name_parts) != 3) or
            len(file_name_parts) == 3 and (file_name_parts[1] != "pdb" or
                                           file_name_parts[2] != "gz")) or \
                (len(file_name_parts) == 2 and file_name_parts[1] != "pdb"):
            print(error_msg + wrong_naming_msg)
            return 1

        file_name_no_ext = file_name_parts[0]
        file_name_no_ext_parts = file_name_no_ext.split(sep="_")
        structure_name = file_name_no_ext_parts[0]
        chain_1 = file_name_no_ext_parts[1]
        chain_2 = file_name_no_ext_parts[2]

        if not structure_name.isalnum():
            print(error_msg + "File name must be alphanumeric")
            return 1

        if not is_chain_in_structure(os.path.join(path, file_name), chain_1):
            print(error_msg + f"Chain {chain_1} not present in the structure")
            return 1

        if not is_chain_in_structure(os.path.join(path, file_name), chain_2):
            print(error_msg + f"Chain {chain_2} not present in the structure")
            return 1

    return 0


def is_chain_in_structure(file_path, chain):
    pdb_parser = PDBParser()
    if file_path.split(sep=".")[-1] == "gz":
        with gzip.open(file_path, 'rt') as pdb_file:
            structure = pdb_parser.get_structure("", pdb_file)
    else:
        structure = pdb_parser.get_structure("", file_path)

    chains = [c.get_id() for c in structure.get_chains()]
    return chain in chains


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-directory", required=True,
                        help="directory containing the input structure files")

    args = parser.parse_args()
    parse_res = parse_input_directory(args.input_directory)
    if parse_res == 1:
        sys.exit(1)
