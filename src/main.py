import argparse
import os
import gzip
import sys
import logging as log
from pdb_processing import *

pdbs = dict()


def parse_input_directory(path):
    # Check if input directory is a directory
    input_dir_error_msg = "Error reading input directory. "
    if not os.path.isdir(path):
        log.error(input_dir_error_msg + f"'{path}' is not a directory")
        sys.exit(1)

    file_names = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
                  and not f.startswith('.')]

    # Check each file
    for file_name in file_names:
        file_name_parts = file_name.split(".")

        # Error messages
        file_error_msg = f"Error in file '{file_name}': "
        wrong_naming_msg = "wrong file naming structure. It should be as follows: " \
                           "<name>_<chain1>_<chain2>.pdb(.gz) "

        # Check extensions and name
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
            log.error(input_dir_error_msg + file_error_msg + ": PDB name must be alphanumeric")
            sys.exit(1)

        # Process PDB
        pdbs[file_name_no_ext] = get_pdb_structure(os.path.join(path, file_name), file_name_no_ext)
        structure = pdbs[file_name_no_ext]

        # Check that the chains in the file name are in the PDB
        chains = [c.get_id() for c in structure.get_chains()]

        if chain_1 not in chains:
            log.error(input_dir_error_msg + file_error_msg + f": chain {chain_1} not present in "
                                                             "the structure")
            sys.exit(1)

        if chain_2 not in chains:
            log.error(input_dir_error_msg + file_error_msg + f": chain {chain_2} not present in "
                                                             "the structure")
            sys.exit(1)

        if len(chains) > 2:
            log.error(input_dir_error_msg + file_error_msg + f": there are more than 2 chains in "
                                                             "the structure")
            sys.exit(1)

    return 0


def get_pdb_structure(file_path, pdb_id):
    pdb_parser = PDBParser()

    if file_path.split(sep=".")[-1] == "gz":
        with gzip.open(file_path, 'rt') as pdb_file:
            return pdb_parser.get_structure(pdb_id, pdb_file)
    else:
        return pdb_parser.get_structure(pdb_id, file_path)


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
    parser.add_argument("-it", "--identity-threshold", default=0.95,
                        help="Minimum percentage of sequence similarity (between 0 and 1) "
                             "to consider two PDB chains the same")
    parser.add_argument("-Rt", "--RMSD-threshold", default=2,
                        help="Maximum RMSD value to consider two (similar) PDB chains the same")
    args = parser.parse_args()

    if args.verbose:
        log.basicConfig(format="%(message)s", level=log.INFO)
    else:
        log.basicConfig(format="%(message)s", level=log.WARNING)

    parse_input_directory(args.input_directory)
    parse_output_directory(args.output_directory, args.force)

    process_pdbs(pdbs, args.identity_threshold, args.RMSD_threshold)


if __name__ == "__main__":
    main()
