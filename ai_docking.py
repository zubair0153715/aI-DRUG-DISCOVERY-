import os
import pickle
import random
from rdkit import Chem
from rdkit.Chem import Descriptors

print("="*60)
print("AI + DOCKING CLOSED-LOOP SYSTEM...")
print("="*60)

# Model load karo
model_file = os.path.join(os.getcwd(), "models", "simple_ai_model.pkl")
if os.path.exists(model_file):
    with open(model_file, "rb") as f:
        model = pickle.load(f)
    fragments = model["fragments"]
    print("✅ AI Model Loaded!")
else:
    print("❌ Model not found. Run train_simple_model.py first.")
    exit()

# Docking Score Simulation (Real docking mein AutoDock Vina use hota hai)
# Yahan hum AI ko batayenge ki kaunsa molecule protein ke saath better bind karega
def calculate_binding_score(smiles):
    """Simple binding score calculator (Real system mein AutoDock Vina use hoga)"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rings = Descriptors.RingCount(mol)
        
        # Better binding ke liye ideal properties
        # (Real docking mein protein-ligand interaction calculate hota hai)
        score = 0
        if 300 < mw < 450: score += 2  # Ideal size
        if 0 < logp < 3: score += 2    # Ideal solubility
        if hbd <= 3: score += 1        # H-bond donors
        if hba <= 6: score += 1        # H-bond acceptors
        if rings >= 1: score += 1      # Ring structures bind better
        
        return score
    return 0

print("\n" + "="*60)
print("AI GENERATING & SCREENING MOLECULES FOR TARGET...")
print("="*60)

best_molecules = []
for i in range(10):
    best_score = -1
    best_smiles = None
    
    # AI 20 molecules generate karta hai
    for attempt in range(20):
        num_frags = random.randint(3, 6)
        smiles = ""
        for _ in range(num_frags):
            smiles += random.choice(fragments)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            score = calculate_binding_score(smiles)
            if score > best_score:
                best_score = score
                best_smiles = smiles
    
    if best_smiles:
        best_molecules.append((best_smiles, best_score))
        print(f"\n🎯 Round {i+1}: Best Molecule Found!")
        print(f"   SMILES: {best_smiles}")
        print(f"   Binding Score: {best_score}/7")

print("\n" + "="*60)
print("TOP 5 DRUG CANDIDATES FOR TARGET PROTEIN:")
print("="*60)

sorted_molecules = sorted(best_molecules, key=lambda x: x[1], reverse=True)[:5]
for i, (smiles, score) in enumerate(sorted_molecules, 1):
    mol = Chem.MolFromSmiles(smiles)
    print(f"\n🏆 Rank #{i}:")
    print(f"   SMILES: {smiles}")
    print(f"   Binding Score: {score}/7")
    print(f"   MW: {Descriptors.MolWt(mol):.2f}")
    print(f"   LogP: {Descriptors.MolLogP(mol):.2f}")

# Save top molecules
output_file = os.path.join(os.getcwd(), "output", "top_drug_candidates.csv")
with open(output_file, "w") as f:
    f.write("Rank,SMILES,Binding_Score,MW,LogP\n")
    for i, (smiles, score) in enumerate(sorted_molecules, 1):
        mol = Chem.MolFromSmiles(smiles)
        f.write(f"{i},{smiles},{score},{Descriptors.MolWt(mol):.2f},{Descriptors.MolLogP(mol):.2f}\n")

print(f"\n✅ Top candidates saved at: {output_file}")
print("="*60)
print("AI DRUG DISCOVERY PIPELINE COMPLETE!")
print("="*60)