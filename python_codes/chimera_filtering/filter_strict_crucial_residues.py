import os
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import BiopythonWarning
import warnings

warnings.simplefilter('ignore', BiopythonWarning)

# Residues likely involved in H-bonding, not hydrophobic
crucial_residues = {
    "SER", "THR", "ASN", "GLN", "ARG", "HIS", "LYS", "ASP", "GLU", "TYR", "CYS"
}
ligand_resname = "7V7"
target_files = {"5tzo.pdb", "8tfq.pdb", "7u64.pdb", "8va0.pdb"}

input_folder = r"C:\Users\Aina Rif'ah\Downloads\UKM\FOR FYPPP\13 PROTEINS PDB"  # <-- update this
output_folder = "refiltered_for_assam"
os.makedirs(output_folder, exist_ok=True)

class StrictHbondResidueSelect(Select):
    def __init__(self, ligand_atoms, max_distance=3.0):
        self.ligand_atoms = ligand_atoms
        self.max_distance = max_distance
        self.selected_residues = set()

    def accept_residue(self, residue):
        if residue.get_resname() not in crucial_residues:
            return 0
        for atom in residue.get_atoms():
            if atom.element not in {"O", "N"}:
                continue
            if atom.get_name() in {"N", "O", "C", "CA"}:  # skip backbone atoms
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

for filename in os.listdir(input_folder):
    if filename not in target_files:
        continue

    path = os.path.join(input_folder, filename)
    structure = parser.get_structure(filename[:4], path)

    ligand_atoms = [atom for model in structure
                         for chain in model
                         for residue in chain
                         if residue.get_resname() == ligand_resname
                         for atom in residue]

    if not ligand_atoms:
        print(f"[SKIPPED] {filename}: No ligand {ligand_resname}")
        continue

    selector = StrictHbondResidueSelect(ligand_atoms)
    io.set_structure(structure)
    filtered_file = os.path.join(output_folder, f"filtered_{filename}")
    io.save(filtered_file, selector)

    n = len(selector.selected_residues)
    if 3 <= n <= 12:
        print(f"[KEPT] {filename}: {n} strict crucial residues")
    else:
        os.remove(filtered_file)
        print(f"[REMOVED] {filename}: {n} residues (still out of range)")
