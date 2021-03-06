import os
import gzip
from Bio.PDB import PDBParser, NeighborSearch, Selection, Polypeptide


def analyse_pdb(pdb_files_path):
    file_names_with_ext = [f for f in os.listdir(pdb_files_path) if
                           os.path.isfile(os.path.join(pdb_files_path, f))
                           and not f.startswith('.')]

    file_names_no_ext = [file_name.split(".")[0] for file_name in file_names_with_ext]

    # structure -> models -> chains -> residues -> atoms
    pdb_structures = [
        get_pdb_structure(os.path.join(pdb_files_path, pdb_name_with_ext), pdb_name_no_ext) for
        pdb_name_with_ext, pdb_name_no_ext in
        zip(file_names_with_ext, file_names_no_ext)]

    structure = pdb_structures[0]
    models = structure.get_models()
    chains = structure.get_chains()
    residues = structure.get_residues()
    atoms = structure.get_atoms()  # this is a generator; call list(atoms) to use as a list
    atoms_list = Selection.unfold_entities(structure, "A")

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


def get_pdb_structure(file_path, pdb_id):
    pdb_parser = PDBParser()

    if file_path.split(sep=".")[-1] == "gz":
        with gzip.open(file_path, 'rt') as pdb_file:
            return pdb_parser.get_structure(pdb_id, pdb_file)
    else:
        return pdb_parser.get_structure(pdb_id, file_path)


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