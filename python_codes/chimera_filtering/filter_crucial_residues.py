import os
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import BiopythonWarning
import warnings

warnings.simplefilter('ignore', BiopythonWarning)

# Only residues likely to form H-bonds (non-hydrophobic)
crucial_residues = {
    "SER", "THR", "ASN", "GLN", "ARG", "HIS", "LYS", "ASP", "GLU", "TYR", "CYS"
}
ligand_resname = "7V7"
input_folder = r"C:\Users\Aina Rif'ah\Downloads\UKM\FOR FYPPP\13 PROTEINS PDB"  # <-- Update this
output_folder = "filtered_crucial_pdbs"
os.makedirs(output_folder, exist_ok=True)

class CrucialResidueSelect(Select):
    def __init__(self, ligand_atoms, max_distance=3.5):
        self.ligand_atoms = ligand_atoms
        self.max_distance = max_distance
        self.selected_residues = set()

    def accept_residue(self, residue):
        if residue.get_resname() not in crucial_residues:
            return 0
        for atom in residue.get_atoms():
            if atom.element not in {"O", "N"}:
                continue
            for ligand_atom in self.ligand_atoms:
                if atom - ligand_atom <= self.max_distance:
                    self.selected_residues.add(residue)
                    return 1
        return 0

    def accept_atom(self, atom):
        return 1 if atom.get_parent() in self.selected_residues else 0

parser = PDBParser(QUIET=True)
io = PDBIO()

total = kept = skipped = deleted = 0

for filename in os.listdir(input_folder):
    if not filename.endswith(".pdb"):
        continue

    total += 1
    path = os.path.join(input_folder, filename)
    structure = parser.get_structure(filename[:4], path)

    ligand_atoms = [atom for model in structure
                         for chain in model
                         for residue in chain
                         if residue.get_resname() == ligand_resname
                         for atom in residue]

    if not ligand_atoms:
        print(f"[SKIPPED] {filename}: No {ligand_resname} found")
        skipped += 1
        continue

    selector = CrucialResidueSelect(ligand_atoms)
    io.set_structure(structure)
    filtered_file = os.path.join(output_folder, f"filtered_{filename}")
    io.save(filtered_file, selector)

    n = len(selector.selected_residues)
    if 3 <= n <= 12:
        print(f"[KEPT] {filename}: {n} crucial residues")
        kept += 1
    else:
        os.remove(filtered_file)
        print(f"[REMOVED] {filename}: {n} residues (out of range)")
        deleted += 1

print("\n--- SUMMARY ---")
print(f"Total files processed: {total}")
print(f"Skipped (no ligand): {skipped}")
print(f"Filtered and kept (3–12 residues): {kept}")
print(f"Removed (too few/many residues): {deleted}")

