from rdkit import Chem
from rdkit.Chem import Descriptors

print("="*50)
print("AI DRUG GENERATOR TEST")
print("="*50)

# Yeh molecules AI ne banaye (Example ke liye)
test_molecules = [
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen (existing drug)
    "CC(=O)Nc1ccc(O)cc1",          # Paracetamol (existing drug)
    "CN1CCC[C@H]1c1cccnc1",       # Nicotine analog
]

print("\nAI Generated Molecules Test:")
print("-"*50)
for smiles in test_molecules:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        print(f"✅ Valid Molecule: {smiles}")
        print(f"   Molecular Weight: {Descriptors.MolWt(mol):.2f}")
        print(f"   LogP: {Descriptors.MolLogP(mol):.2f}")
        print(f"   H-Bond Donors: {Descriptors.NumHDonors(mol)}")
        print(f"   H-Bond Acceptors: {Descriptors.NumHAcceptors(mol)}")
        print()
    else:
        print(f"❌ Invalid: {smiles}")

print("="*50)
print("AI SYSTEM READY FOR DRUG GENERATION!")
print("="*50)