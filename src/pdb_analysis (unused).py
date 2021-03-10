import os
import gzip
from Bio.PDB import PDBParser, NeighborSearch, Selection, Polypeptide
from Bio.pairwise2 import align, format_alignment


def analyse_pdb(pdb_files_path):
    # structure -> models -> chains -> residues -> atoms
    pdb_structures = get_pdb_structures(pdb_files_path)
    pdb_ids = list(pdb_structures.keys())

    structure = pdb_structures[pdb_ids[0]]
    models = structure.get_models()
    chains = structure.get_chains()
    residues = structure.get_residues()
    atoms = structure.get_atoms()  # this is a generator; call list(atoms) to use as a list
    atoms_list = Selection.unfold_entities(structure, "A")

    struct_A_B = pdb_structures[pdb_ids[0]]
    struct_A_C = pdb_structures[pdb_ids[1]]

    chains_A_B = struct_A_B.get_chains()
    chains_A_C = struct_A_C.get_chains()

    struct_from_chain = Selection.unfold_entities(next(chains_A_B), "S")

    model = structure[0]
    chain = model["A"]

    residue = chain[100]
    residue_name = residue.get_resname()
    has_atom_ca = residue.has_id("CA")
    print("Is the residue an amino acid? ", Polypeptide.is_aa(residue))

    atom = residue["CA"]
    atom_coords = atom.get_coord()
    atom_id = atom.get_id()

    residue2 = chain[101]
    atom2 = residue2["CA"]
    distance = atom - atom2
    print(f"Distance between atoms {atom_id} and {atom2.get_id()} = {distance} Å")

    neighbor_search = NeighborSearch(list(atoms))
    radius = 2
    neighbor_atoms = neighbor_search.search(atom_coords, radius)
    print(f"{len(neighbor_atoms) - 1} atoms are in a radius of {radius} Å")

    ppb = Polypeptide.PPBuilder()
    ppb_polypeptides = ppb.build_peptides(structure)
    ppb_sequences = [pp.get_sequence() for pp in ppb_polypeptides]

    cappb = Polypeptide.CaPPBuilder()
    cappb_polypeptides = cappb.build_peptides(structure)
    sequences = [pp.get_sequence() for pp in cappb_polypeptides]

    print("Are C-N and CA-CA polypeptides builds equals? ", ppb_sequences == sequences)

    alignments = align.globalxx(sequences[0], sequences[1], one_alignment_only=True)
    # print(format_alignment(*alignments[0]))

    sequences_by_chain = dict(
        (struct_id, get_sequences_by_chain(struct)) for struct_id, struct in pdb_structures.items())


def get_pdb_structures(pdb_files_path):
    file_names_with_ext = [f for f in os.listdir(pdb_files_path) if
                           os.path.isfile(os.path.join(pdb_files_path, f))
                           and not f.startswith('.')]

    file_names_no_ext = [file_name.split(".")[0] for file_name in file_names_with_ext]

    return dict((pdb_name_no_ext,
                 get_pdb_structure(os.path.join(pdb_files_path, pdb_name_with_ext),
                                   pdb_name_no_ext)) for
                pdb_name_with_ext, pdb_name_no_ext in
                zip(file_names_with_ext, file_names_no_ext))


def get_pdb_structure(file_path, pdb_id):
    pdb_parser = PDBParser()

    if file_path.split(sep=".")[-1] == "gz":
        with gzip.open(file_path, 'rt') as pdb_file:
            return pdb_parser.get_structure(pdb_id, pdb_file)
    else:
        return pdb_parser.get_structure(pdb_id, file_path)


def get_sequences_by_chain(structure):
    cappb = Polypeptide.CaPPBuilder()
    cappb_polypeptides = cappb.build_peptides(structure)

    return dict((chain.get_id(), pp.get_sequence()) for pp, chain in
                zip(cappb_polypeptides, structure.get_chains()))


def compare_sequences(seq1, seq2, threshold=0.95):
    alignment = align.globalxx(seq1, seq2, one_alignment_only=True)[0]
    if (alignment.score / len(seq1)) > threshold:  # or max(len(seq1), len(seq2))?
        return True
    return False


def print_structure(structure):
    print("Structure ", structure.get_id())
    for model in structure:
        print("Model", model.get_id())
        for chain in model:
            print("Chain", chain.get_id())
            for residue in chain:
                print("Residue", residue)
                for atom in residue:
                    print("Atom", atom)
