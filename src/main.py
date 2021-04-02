import argparse
import os
import gzip
import sys
import logging as log
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1
from pdb_processing import process_pdbs
from reconstruction import build_complex

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
        # Error messages
        file_error_msg = f"Error in file '{file_name}': "
        wrong_naming_msg = "wrong file naming structure. It should be as follows: " \
                           "<PDB name>_<chain1>_<chain2>.pdb(.gz) or " \
                           "<Uniprot ID>.<DNA or RNA>.<PDB name>_<chain1>_<chain2>.pdb(.gz)"

        file_name_parts = file_name.split(".")

        # Check extensions and name
        has_prot_prot_format_len = len(file_name_parts) == 2 or len(file_name_parts) == 3
        has_prot_dna_format_len = len(file_name_parts) == 4 or len(file_name_parts) == 5

        has_pdb_format = file_name_parts[-1] == "pdb"
        has_compressed_format = file_name_parts[-1] == "gz" and file_name_parts[-2] == "pdb"

        if not (has_prot_prot_format_len or has_prot_dna_format_len) or (
                len(file_name_parts) == 2 and not has_pdb_format) or (
                len(file_name_parts) == 3 and not has_compressed_format) or (
                len(file_name_parts) == 4 and not has_pdb_format) or (
                len(file_name_parts) == 5 and not has_compressed_format):
            log.error(input_dir_error_msg + file_error_msg + wrong_naming_msg)
            sys.exit(1)

        file_name_parts_no_ext = [part for part in file_name_parts
                                  if part != "pdb" and part != "gz"]

        structure_name = file_name_parts_no_ext[-1]
        structure_name_parts = structure_name.split(sep="_")
        pdb_name = structure_name_parts[0]
        chains_in_file_name = structure_name_parts[1:]

        # Only DNA has double chain
        '''
        if has_prot_dna_format_len and file_name_parts_no_ext[1] == "DNA":
            split_double_chain = list(chains_in_file_name[1])
            chains_in_file_name[1] = split_double_chain[0]
            chains_in_file_name.append(split_double_chain[1])
        '''

        if has_prot_dna_format_len:  # Don't check if [1] is DNA because it can be RNA??
            chains_in_file_name = [chain_letter for chain_letter in "".join(chains_in_file_name)]

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
            if chain not in chains:
                log.error(input_dir_error_msg + file_error_msg + f": chain {chain} not present"
                                                                 " in the structure")

    return 0


def get_pdb_structure(file_path, pdb_id):
    if file_path.split(sep=".")[-1] == "gz":
        with gzip.open(file_path, 'rt') as pdb_file:
            structure = PDBParser().get_structure(pdb_id, pdb_file)
            add_chain_sequences(structure)
            return structure
    else:
        structure = PDBParser().get_structure(pdb_id, file_path)
        add_chain_sequences(structure)
        return structure


def add_chain_sequences(structure):
    for chain in structure.get_chains():
        # Look at the first non-heteroatom residue to see if it is RNA, DNA or a protein
        first_residue = next(res for res in chain.get_residues() if res.get_id()[0].isspace())
        if first_residue.get_resname().startswith(" "):
            # DNA or RNA sequence with 1 letter ([-1] position of resname)
            chain.xtra["type"] = "nuc"
            chain_sequence = "".join([res.get_resname()[-1] for res in chain.get_residues()
                                      if res.get_id()[0].isspace()])
        else:
            # Protein sequence
            chain.xtra["type"] = "prot"
            chain_sequence = "".join([protein_letters_3to1[res.get_resname().lower().capitalize()]
                                      for res in chain.get_residues() if res.get_id()[0].isspace()])
        chain.xtra["seq"] = chain_sequence


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
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser.add_argument("-Rt", "--RMSD-threshold", default=2.5,
                        help="Maximum RMSD value to consider two (similar) PDB chains the same")
    parser.add_argument("-ns", "--Neighbor-Search-distance", default=3.5,
                        help="Minimum distance between two PDB chains to consider that "
                             "they are actually interacting")
    parser.add_argument("-cd", "--clashes-distance", default=1,
                        help="Maximum distance between two PDB interacting chains to consider that "
                             "they have clashes between them")
    parser.add_argument("-cad", "--CA-atoms-distance", default=5,
                        help="Maximum distance between Alpha carbons atoms of two PDB interacting "
                             "chains to consider they have clashes between them")
    parser.add_argument("-nc", "--number-clashes", default=10,
                        help="Maximum number of close atoms to consider two PDB interacting "
                             "chains can have a clash")
    args = parser.parse_args()

    if args.verbose:
        log.basicConfig(format="%(message)s", level=log.INFO)
    else:
        log.basicConfig(format="%(message)s", level=log.ERROR)
        log.captureWarnings(True)

    parse_input_directory(args.input_directory)
    parse_output_directory(args.output_directory, args.force)

    process_pdbs(pdbs, args.identity_threshold, args.Neighbor_Search_distance, args.RMSD_threshold)

    build_complex(args.output_directory, args.clashes_distance, args.CA_atoms_distance,
                  args.number_clashes)


if __name__ == "__main__":
    main()
