import argparse
import os
import gzip
import sys
import logging as log
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1 as Res_dict
from PDB_processing import process_pdbs
from Construction import build_complex
from PDB_tools import get_chain_full_id, minimize
from Stoichiometry import *

pdb_chains = list()
stoich_dict = dict()


def parse_input_directory(path, stoichiometry_path):
    log.info("Parsing input PDB files")
    # Check if input directory is a directory
    input_dir_error_msg = "Error reading input directory. "
    if not os.path.isdir(path):
        log.error(input_dir_error_msg + f"'{path}' is not a directory")
        sys.exit(1)

    file_names = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
                  and not f.startswith('.')]

    all_chains = set()
    all_chains_by_pdb = dict()
    all_prefixes = set()

    # Check each file
    for file_name in file_names:
        # Error messages
        file_error_msg = f"Error in file '{file_name}': "
        wrong_naming_msg = "wrong file naming structure. It should be as follows\n " \
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
        if len(structure_name_parts) != 3:
            log.error(input_dir_error_msg + file_error_msg + wrong_naming_msg)
            sys.exit(1)
        pdb_name = structure_name_parts[0]
        chains_in_file_name = structure_name_parts[1:]

        if has_prot_dna_format_len:
            chains_in_file_name = [chain_letter for chain_letter in "".join(chains_in_file_name)]

        if len(pdb_name) != 4:
            log.error(input_dir_error_msg + file_error_msg + "PDB name must have 4 characters:\n"
                                "<PDB name>_<chain1>_<chain2>.pdb(.gz) or "
                                "<Uniprot ID>.<DNA or RNA>.<PDB name>_<chain1>_<chain2>.pdb(.gz)")
            sys.exit(1)

        if not pdb_name.isalnum():
            log.error(input_dir_error_msg + file_error_msg + "PDB name must be alphanumeric")
            sys.exit(1)

        # Check prefix
        if len(file_name_parts_no_ext) > 1:
            prefix = file_name_parts_no_ext[0]
            if len(prefix) != 6 or not prefix[0].isalnum() or not prefix[1:].isnumeric():
                log.error(input_dir_error_msg + file_error_msg +
                          "Uniprot ID must have one letter and five numbers")
                sys.exit(1)
            all_prefixes.add(prefix)

        # Process PDB
        file_name_no_ext = ".".join(file_name_parts_no_ext)
        structure = get_pdb_structure(os.path.join(path, file_name), file_name_no_ext)

        # Add chains to the pdb_chains list
        for chain in structure.get_chains():
            pdb_chains.append(chain)

        if pdb_name not in all_chains_by_pdb:
            all_chains_by_pdb[pdb_name] = set()

        # Check that the chains in the file name are in the PDB
        chain_ids = [c.get_id() for c in structure.get_chains()]
        for chain_id in chains_in_file_name:
            if chain_id not in chain_ids:
                log.error(input_dir_error_msg + file_error_msg + f"chain {chain_id} not present"
                                                                 " in the structure")
                sys.exit(1)
            else:
                all_chains.add(chain_id)
                all_chains_by_pdb[pdb_name].add(chain_id)
                chain_ids.remove(chain_id)


    if stoichiometry_path:
        parse_stoichiometry(stoichiometry_path, all_chains_by_pdb, all_prefixes, all_chains, stoich_dict)

    return 0



def get_pdb_structure(file_path, pdb_id):
    if file_path.split(sep=".")[-1] == "gz":
        with gzip.open(file_path, 'rt') as pdb_file:
            structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)
            add_chain_sequences(structure)
            return structure
    else:
        structure = PDBParser(QUIET=True).get_structure(pdb_id, file_path)
        add_chain_sequences(structure)
        return structure



def add_chain_sequences(structure):
    for chain in structure.get_chains():
        # Look at the first non-heteroatom residue to see if it is RNA, DNA or a protein
        first_residue = next(res for res in chain.get_residues() if res.get_id()[0].isspace())
        if first_residue.get_resname().startswith(" "):
            # DNA or RNA sequence with 1 letter ([-1] position of resname)
            chain.xtra["type"] = "nuc"
            # To check if the residue is a ribonucleotide ("  X") or a deoxyribonucleotide (" DX")
            nuc_letters = "A G T C U".split(" ")
            is_ribonuc = lambda res: res.get_resname()[:-1] == "  " and \
                                     res.get_resname()[-1] in nuc_letters
            is_deoxyribonuc = lambda res: res.get_resname()[:-1] == " D" and \
                                          res.get_resname()[-1] in nuc_letters
            # Get the sequence
            chain_sequence = "".join([res.get_resname()[-1] for res in chain.get_residues()
                                      if (is_ribonuc(res) or is_deoxyribonuc(res))])
        else:
            # Protein sequence
            chain.xtra["type"] = "prot"
            # Function to change the PDB names of the residues from RES (e.g., HIS) to Res (His)
            # to query the residue dictionary (keys as "Res")
            RES_to_Res = lambda res: res.get_resname().lower().capitalize()
            # Get the sequence
            chain_sequence = "".join([Res_dict[RES_to_Res(res)] for res in chain.get_residues()
                                      if RES_to_Res(res) in Res_dict])
        chain.xtra["seq"] = chain_sequence
        # Also save the original chain ID as <PDB file name>:<Chain ID>
        chain.xtra["full_id"] = get_chain_full_id(chain)



def parse_output_directory(path_dir, force):
    output_dir_error_msg = "Error with output directory. "
    if force is True:
        log.info("Removing files from output directories")
        for folder in ['structures', 'analysis']:
            dir = os.path.join(path_dir, folder)
            os.makedirs(dir, exist_ok=True)
            for file in os.listdir(dir):
                os.remove(os.path.join(dir, file))
    else:
        if not os.path.exists(path_dir):
            log.info("Creating output directories")
            for folder in ['structures', 'analysis']:
                os.makedirs(os.path.join(path_dir, folder), exist_ok=False)
        else:
            log.error(output_dir_error_msg + f"'{path_dir}' directory already exists. Use [-f, "
                                             "--force] option to overwrite.")
            sys.exit(1)
    return 0


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""MacroPy 1.0 -
                                                 Reconstruct a whole biological macro-complex using PDBs of
                                                 its pairwise interactions as input, either protein-protein,
                                                 protein-DNA or protein-RNA.""")
    parser.add_argument("-i", "--input-directory", required=True,
                        help="Directory containing the input structure files")
    parser.add_argument("-o", "--output-directory", required=True,
                        help="Create the output directories")
    parser.add_argument("-c", "--complex-name", default="Complex",
                        help="Reconstructed complex name")
    parser.add_argument("-f", "--force", action="store_true", default=False,
                        help="Force overwriting if the output directory already exists")
    parser.add_argument("-s", "--stoichiometry",
                        help="File containing the stoichiometry")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="Program log will be printed to standard error while running")
    parser.add_argument("-min", "--minimization", nargs='?', default=False, const=True,
                        help="Perform an energy minimization by Conjugate Gradients algorithm with the specified "
                             "(-min X) number of steps (or, if no number is specified (-min), until convergence)")
    parser.add_argument("-pdb", "--save-pdb", action="store_true", default=False,
                        help="Besides the .cif file, save a .pdb file with up to 25 chains")
    parser.add_argument("-mc", "--max-chains", default=180,
                        help="Number of chains of the complex at which to stop adding new chains")
    parser.add_argument("-it", "--identity-threshold", default=0.95,
                        help="Minimum percentage of sequence similarity (between 0 and 1) "
                             "to consider two PDB chains the same")
    parser.add_argument("-Rt", "--RMSD-threshold", default=2.5,
                        help="Maximum RMSD value to consider two (similar) PDB chains the same")
    parser.add_argument("-ns", "--Neighbor-Search-distance", default=5,
                        help="Minimum distance between two PDB chains to consider that "
                             "they are actually interacting")
    parser.add_argument("-cd", "--clashes-distance", default=1.8,
                        help="Maximum distance between atoms of two chains to consider that "
                             "they have clashes between them")
    parser.add_argument("-nc", "--number-clashes", default=20,
                        help="Maximum number of close atoms to consider that two chains "
                             " are clashing")
    args = parser.parse_args()

    if args.verbose:
        log.basicConfig(format="%(message)s", level=log.INFO)
    else:
        log.basicConfig(format="%(message)s", level=log.ERROR)
        log.captureWarnings(True)

    log.info("Running MacroPy")
    parse_input_directory(args.input_directory, args.stoichiometry)
    parse_output_directory(args.output_directory, args.force)

    process_pdbs(pdb_chains, args.identity_threshold, args.Neighbor_Search_distance, args.RMSD_threshold)

    build_complex(args.output_directory, int(args.max_chains), args.number_clashes, args.save_pdb,
                  args.clashes_distance, stoich_dict, str(args.complex_name), args.identity_threshold)

    if args.minimization:
        minimize(args.output_directory, args.complex_name, args.minimization)


if __name__ == "__main__":
    main()
