from Bio.PDB import PDBParser, PDBIO, Select

# Input/output file names
input_pdb = "8ef5_structure.pdb"  # Change this to every original PDB filename
output_pdb = "filtered_12_residues.pdb"

ligand_resname = "7V7"
distance_cutoff = 4.0

# Hydrophobic residues to exclude
hydrophobic_residues = ["ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "MET", "TRP"]

parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", input_pdb)

# Get ligand atoms
ligand_atoms = []
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_resname() == ligand_resname:
                ligand_atoms.extend(residue.get_atoms())

# Find residues within distance cutoff and exclude hydrophobic
selected_residues = []
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_resname() == ligand_resname:
                continue  # skip ligand itself
            # check distance to any ligand atom
            close = False
            for atom in residue:
                for ligand_atom in ligand_atoms:
                    if (atom - ligand_atom) < distance_cutoff:
                        close = True
                        break
                if close:
                    break
            if close and residue.get_resname() not in hydrophobic_residues:
                selected_residues.append(residue)

# Limit to max 12 residues
selected_residues = selected_residues[:12]

class ResidueSelect(Select):
    def accept_residue(self, residue):
        return residue in selected_residues

io = PDBIO()
io.set_structure(structure)
io.save(output_pdb, ResidueSelect())

print(f"Filtered structure with max 12 residues saved to {output_pdb}")
