import argparse
import os
import gzip
import sys
import logging as log
from pdb_processing import *
from biskit import PDBModel

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
                           "<PDB name>_<chain1>_<chain2>.pdb(.gz) or " \
                           "<Uniprot ID>.<DNA or RNA>.<PDB name>_<chain1>_<chain2>.pdb(.gz)"

        # Check extensions and name
        parts = file_name_parts
        if (len(parts) > 3 and \
            (( parts[-1] == "gz" and (len(parts) != 5 or parts[-2] != "pdb") ) or \
            ( parts[-1] == "pdb" and len(parts) != 4 ))) or \
           (len(parts) <= 3 and \
            (( parts[-1] == "gz" and (len(parts) != 3 or parts[-2] != "pdb") ) or \
             ( parts[-1] == "pdb" and len(parts) != 2 ))):

            # if (len(file_name_parts) != 2 and len(file_name_parts) != 3) or \
            #         (len(file_name_parts) == 3 and
            #          (file_name_parts[1] != "pdb" or file_name_parts[2] != "gz")) or \
            #         (len(file_name_parts) == 2 and file_name_parts[1] != "pdb"):
            log.error(input_dir_error_msg + file_error_msg + wrong_naming_msg)
            sys.exit(1)

        file_name_parts_no_ext = [ part for part in file_name_parts if part != "pdb" and part != "gz" ]

        structure_name = file_name_parts_no_ext[-1]
        structure_name_parts = structure_name.split(sep="_")
        pdb_name = structure_name_parts.pop(0)
        chains_in_file_name = structure_name_parts

        if not pdb_name.isalnum():
            log.error(input_dir_error_msg + file_error_msg + ": PDB name must be alphanumeric")
            sys.exit(1)

        # Process PDB
        file_name_no_ext = ".".join(file_name_parts_no_ext)
        pdbs[file_name_no_ext] = get_pdb_structure(os.path.join(path, file_name), file_name_no_ext)
        structure = pdbs[file_name_no_ext]

        # Check that the chains in the file name are in the PDB
        chains = [c.get_id() for c in structure.get_chains()]

        for chain in chains_in_file_name:
            for letter in chain:
                if letter not in chains:
                    log.error(input_dir_error_msg + file_error_msg + f": chain {letter} not present"
                                                             " in the structure")
                    sys.exit(1)
    return 0



def get_pdb_structure(file_path, pdb_id):

    if file_path.split(sep=".")[-1] == "gz":
        with gzip.open(file_path, 'rt') as pdb_file:
            structure = PDBParser().get_structure(pdb_id, pdb_file)
            structure.xtra = get_biskit_seqs(structure, pdb_file, pdb_id)
            return structure
    else:
        structure = PDBParser().get_structure(pdb_id, file_path)
        structure.xtra = get_biskit_seqs(structure, file_path, pdb_id)
        return structure



def get_biskit_seqs(structure, pdb_file, pdb_id):
    chain_list = list(structure.get_chains())
    chain_ids = list(chain.get_id() for chain in chain_list)

    biskit_obj = PDBModel(pdb_file, pdb_id)
    whole_sequence = biskit_obj.sequence()

    atom_chain_breaks = biskit_obj.chainEndIndex()
    res_chain_breaks = list( biskit_obj.atom2resIndices(atom_chain_breaks) )

    chain_seqs = dict()
    prev = 0
    for i in range(0,len(res_chain_breaks)):
        ix = list(res_chain_breaks + [None])[i]
        chain_seqs[chain_ids[i]] = whole_sequence[prev:ix+1]
        prev = ix + 1

    return chain_seqs




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
