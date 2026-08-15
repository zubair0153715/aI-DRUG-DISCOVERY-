import os
import pickle
import random
from rdkit import Chem
from rdkit.Chem import Descriptors

print("="*60)
print("AI MODEL SE MOLECULES GENERATE KAR RAHE HAIN...")
print("="*60)

# Model load karo
model_file = os.path.join(os.getcwd(), "models", "simple_ai_model.pkl")
with open(model_file, "rb") as f:
    model = pickle.load(f)

fragments = model["fragments"]
print(f"✅ Model Loaded! (Version: {model['version']})")
print(f"✅ Fragments available: {len(fragments)}")
print(f"✅ Previously learned molecules: {len(model['trained_molecules'])}")

print("\n" + "="*60)
print("GENERATING 20 NOVEL DRUG-LIKE MOLECULES...")
print("="*60)

valid_count = 0
generated = []

for i in range(20):
    for attempt in range(50):
        num_frags = random.randint(3, 6)
        smiles = ""
        for _ in range(num_frags):
            smiles += random.choice(fragments)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            if 100 < mw < 500 and -2 < logp < 5 and hbd <= 5 and hba <= 10:
                valid_count += 1
                generated.append(smiles)
                print(f"\n💊 AI Molecule #{valid_count}: {smiles}")
                print(f"   MW: {mw:.2f} g/mol")
                print(f"   LogP: {logp:.2f}")
                print(f"   HBD: {hbd} | HBA: {hba}")
                print(f"   ✅ Drug-Like: YES")
                break

# Generated molecules save karo
output_file = os.path.join(os.getcwd(), "output", "generated_molecules.csv")
with open(output_file, "w") as f:
    f.write("SMILES,MW,LogP,HBD,HBA\n")
    for smiles in generated:
        mol = Chem.MolFromSmiles(smiles)
        f.write(f"{smiles},{Descriptors.MolWt(mol):.2f},{Descriptors.MolLogP(mol):.2f},{Descriptors.NumHDonors(mol)},{Descriptors.NumHAcceptors(mol)}\n")

print("\n" + "="*60)
print(f"RESULT: {valid_count}/20 Valid Drug-Like Molecules Generated!")
print(f"✅ Saved at: {output_file}")
print("AI MODEL FULLY OPERATIONAL!")
print("="*60)