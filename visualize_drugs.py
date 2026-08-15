from rdkit import Chem
from rdkit.Chem import Draw
import os

print("="*50)
print("DRAWING AI GENERATED MOLECULES...")
print("="*50)

# AI ke banaye kuch molecules
molecules_smiles = [
    ("Drug 1", "CC(=O)Nc1ccc(O)cc1"),
    ("Drug 2", "O=C(Nc1ccccc1)c1ccco1"),
    ("Drug 3", "CN1CCN(CC1)c1ccc(cc1)F"),
]

mols = []
legends = []

for name, smiles in molecules_smiles:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mols.append(mol)
        legends.append(name)

# 2D structure draw karo
img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(300,300), legends=legends)

# Save karo
output_path = os.path.join(os.getcwd(), "output", "ai_drugs.png")
img.save(output_path)

print(f"✅ Images saved at: {output_path}")
print("="*50)