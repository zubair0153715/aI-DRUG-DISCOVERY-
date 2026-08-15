import os
import random
from rdkit import Chem
from rdkit.Chem import Descriptors
import pickle

print("="*60)
print("TRAINING SIMPLE AI DRUG GENERATOR...")
print("="*60)

# Common drug fragments (AI inhi se seekhega)
fragments = [
    "c1ccccc1",      # Benzene ring
    "c1ccncc1",      # Pyridine ring
    "c1ccccc1O",     # Phenol
    "c1ccccc1F",     # Fluorobenzene
    "c1ccccc1Cl",    # Chlorobenzene
    "C(=O)N",        # Amide
    "C(=O)O",        # Acid
    "S(=O)(=O)N",    # Sulfonamide
    "C#N",           # Nitrile
    "C(F)(F)F",      # Trifluoromethyl
    "OC",            # Hydroxyl
    "COC",           # Methoxy
    "N",             # Amine
    "CC(=O)",        # Acetyl
    "CN(C)",         # Dimethylamine
    "CC(C)",         # Isopropyl
    "CC",            # Ethyl
    "C",             # Methyl
    "O",             # Oxygen linker
    "NCC",           # Ethylamine
]

print(f"✅ Training Data: {len(fragments)} chemical fragments loaded!")
print("Generating and validating molecules...")

# Model ka "brain" - yeh sikh raha hai kaise molecules banane hain
trained_knowledge = []
attempts = 0
while len(trained_knowledge) < 50 and attempts < 1000:
    attempts += 1
    # Random molecule generate karo
    num_frags = random.randint(3, 6)
    smiles = ""
    for _ in range(num_frags):
        smiles += random.choice(fragments)
    
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        if 100 < mw < 500 and -2 < logp < 5:
            trained_knowledge.append(smiles)

# Model save karo
model_file = os.path.join(os.getcwd(), "models", "simple_ai_model.pkl")
with open(model_file, "wb") as f:
    pickle.dump({
        "fragments": fragments,
        "trained_molecules": trained_knowledge,
        "version": "1.0"
    }, f)

print(f"\n✅ AI Model Trained!")
print(f"✅ Valid molecules learned: {len(trained_knowledge)}")
print(f"✅ Model saved at: {model_file}")
print("="*60)