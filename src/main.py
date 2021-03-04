import argparse
import os
import gzip
import sys
import logging as log
from Bio.PDB import PDBParser
from pdb_analysis import analyse_pdb


def parse_input_directory(path):
    input_dir_error_msg = "Error reading input directory. "
    if not os.path.isdir(path):
        log.error(input_dir_error_msg + f"'{path}' is not a directory")
        sys.exit(1)

    file_names = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
                  and not f.startswith('.')]
    for file_name in file_names:
        file_name_parts = file_name.split(".")
        file_error_msg = f"Error in file '{file_name}': "
        wrong_naming_msg = "wrong file naming structure. It should be as follows: " \
                           "<name>_<chain1>_<chain2>.pdb(.gz) "

        if (len(file_name_parts) != 2 and len(file_name_parts) != 3) or \
                (len(file_name_parts) == 3 and
                 (file_name_parts[1] != "pdb" or file_name_parts[2] != "gz")) or \
                (len(file_name_parts) == 2 and file_name_parts[1] != "pdb"):
            log.error(input_dir_error_msg + file_error_msg + wrong_naming_msg)
            sys.exit(1)

        file_name_no_ext = file_name_parts[0]
        file_name_no_ext_parts = file_name_no_ext.split(sep="_")
        structure_name = file_name_no_ext_parts[0]
        chain_1 = file_name_no_ext_parts[1]
        chain_2 = file_name_no_ext_parts[2]

        if not structure_name.isalnum():
            log.error(input_dir_error_msg + file_error_msg + ": file name must be alphanumeric")
            sys.exit(1)

        if not is_chain_in_structure(os.path.join(path, file_name), chain_1):
            log.error(input_dir_error_msg + file_error_msg + f": chain {chain_1} not present in "
                                                             "the structure")
            sys.exit(1)

        if not is_chain_in_structure(os.path.join(path, file_name), chain_2):
            log.error(input_dir_error_msg + file_error_msg + f": chain {chain_2} not present in "
                                                             "the structure")
            sys.exit(1)

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


def parse_output_directory(path_dir, force):
    output_dir_error_msg = "Error with output directory. "
    if force is True:
        for folder in ['structures', 'analysis']:
            os.makedirs(os.path.join(path_dir, folder), exist_ok=True)
    else:
        if not os.path.exists(path_dir):
            for folder in ['structures', 'analysis']:
                os.makedirs(os.path.join(path_dir, folder), exist_ok=False)
        else:
            log.error(output_dir_error_msg + f"'{path_dir}' directory already exists. Use [-f, "
                                             "--force] option to overwrite.")
            sys.exit(1)
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-directory", required=True,
                        help="Directory containing the input structure files")
    parser.add_argument("-o", "--output-directory", required=True,
                        help="Create the output directories")
    parser.add_argument("-f", "--force", action="store_true", default=False,
                        help="Force overwriting if the output directory already exists")
    parser.add_argument("-s", "--stoichiometry",
                        help="(Directory)/File containing the stoichiometry (file)")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Program log will be printed to standard error while running")
    args = parser.parse_args()

    if args.verbose:
        log.basicConfig(format="%(message)s", level=log.INFO)
    else:
        log.basicConfig(format="%(message)s", level=log.WARNING)

    parse_input_directory(args.input_directory)
    parse_output_directory(args.output_directory, args.force)

    analyse_pdb(args.input_directory)


if __name__ == "__main__":
    main()
