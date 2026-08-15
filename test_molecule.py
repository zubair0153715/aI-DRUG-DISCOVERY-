from rdkit import Chem
from rdkit.Chem import Descriptors

# Yeh Aspirin (Dard ki dawai) ka chemical code hai
smiles = "CC(=O)Oc1ccccc1C(=O)O"

# RDKit is code ko padh kar 3D molecule banayega
mol = Chem.MolFromSmiles(smiles)

if mol is not None:
    print("SUCCESS: Molecule Valid hai!")
    print("Molecular Weight:", Descriptors.MolWt(mol))
    print("LogP (Fat solubility):", Descriptors.MolLogP(mol))
else:
    print("ERROR: Invalid SMILES")