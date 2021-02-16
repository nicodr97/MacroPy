import argparse
import os
from Bio.PDB import PDBParser
import gzip
import sys


def parse_input_directory(path):
    if not os.path.isdir(path):
        print("Path is not a directory") ############ Debería ser un sys.stderr.write mejor?
        return 1 ################### Al final sí que lo cambiamos por sys.exit no?

    file_names = [f for f in os.listdir(path) if
                  os.path.isfile(os.path.join(path, f)) and not f.startswith('.')]
    for file_name in file_names:
        file_name_parts = file_name.lower().rsplit(".", 2) ########## Quizás habría que añadir el lower() para evitar movidas?
        error_msg = f"Error in file '{file_name}'. "
        wrong_naming_msg = "Wrong file naming structure. It should be as follows: <name>_<chain1>_<chain2>.pdb(.gz)"

        if (len(file_name_parts) != 2 and len(file_name_parts) != 3) or \ # There are more than 1 or 2 extensions OR
           \
           (len(file_name_parts) == 3 and \ # There are two extensions 
                (file_name_parts[1] != "pdb" or file_name_parts[2] != "gz")) or \ # but the 1st one isn't pdb or the 2nd one isn't gz OR
           \
           (len(file_name_parts) == 2 and file_name_parts[1] != "pdb"): # There is one extension but it isn't pdb
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


def parse_output_directory(path_dir, force):
    if force is True:
        for i in ['structures','analysis']:
            os.makedirs(os.path.join(path_dir, i), exist_ok=True)
    else:
        if not os.path.exists(path_dir):
            for i in ['structures','analysis']:
                os.makedirs(os.path.join(path_dir, i), exist_ok=False)
        else:
            print(f"{path_dir} directory already exists. Use -f, --force option to overwrite.")
            return 1
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-directory", required=True,
                        help="Directory containing the input structure files")
    parser.add_argument("-o", "--output-directory", required=True,
                        help="Create the output directories")
    parser.add_argument("-f", "--force", action = "store_true", default = False,
                        help="Force overwriting if the ouput directory already exists")
    parser.add_argument("-s", "--stoichiometry",
                        help="(Directory)/File containing the stoichiometry (file)")
    parser.add_argument("-v", "--verbose", action = "store_true", default = False,
                        help="Program log will be printed to standard error while running")
    args = parser.parse_args()

    parse_input_directory(args.input_directory)
    parse_output_directory(args.output_directory, args.force)


    if args.verbose == True:
        outflow = sys.stderr
    else:
        outflow = open(os.devnull, "w")
    
    # Programa

    if outflow != sys.stderr:
        outflow.close()