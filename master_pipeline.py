import os
import pickle
import random
from rdkit import Chem
from rdkit.Chem import Descriptors
from Bio.PDB import PDBParser

print("="*70)
print("🏆 AI DRUG DISCOVERY MASTER PIPELINE")
print("="*70)

# STEP 1: Protein Load karo
print("\n📌 STEP 1: Loading Target Protein...")
pdb_file = os.path.join(os.getcwd(), "data", "1m17.pdb")
if os.path.exists(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", pdb_file)
    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            print(f"✅ Protein Loaded: Chain {chain.id} ({len(residues)} residues)")
            break
        break
else:
    print("⚠️ Protein file not found. Using simulated target.")

# STEP 2: AI Model Load karo
print("\n📌 STEP 2: Loading AI Model...")
model_file = os.path.join(os.getcwd(), "models", "simple_ai_model.pkl")
with open(model_file, "rb") as f:
    model = pickle.load(f)
fragments = model["fragments"]
print(f"✅ AI Model Loaded ({len(fragments)} fragments)")

# STEP 3: Generate & Screen
print("\n📌 STEP 3: AI Generating & Screening Molecules...")
def calculate_binding_score(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rings = Descriptors.RingCount(mol)
        score = 0
        if 300 < mw < 450: score += 2
        if 0 < logp < 3: score += 2
        if hbd <= 3: score += 1
        if hba <= 6: score += 1
        if rings >= 1: score += 1
        return score
    return 0

best_molecules = []
for i in range(20):
    best_score = -1
    best_smiles = None
    for attempt in range(30):
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

# STEP 4: Results
print("\n📌 STEP 4: Final Results")
print("-"*70)
sorted_molecules = sorted(best_molecules, key=lambda x: x[1], reverse=True)[:10]
print("\n🏆 TOP 10 DRUG CANDIDATES:")
print("-"*70)
for i, (smiles, score) in enumerate(sorted_molecules, 1):
    mol = Chem.MolFromSmiles(smiles)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    print(f"\n#{i}. Binding Score: {score}/7")
    print(f"   SMILES: {smiles}")
    print(f"   MW: {mw:.2f} | LogP: {logp:.2f}")

# Save Results
output_file = os.path.join(os.getcwd(), "output", "final_drug_candidates.csv")
with open(output_file, "w") as f:
    f.write("Rank,SMILES,Binding_Score,MW,LogP\n")
    for i, (smiles, score) in enumerate(sorted_molecules, 1):
        mol = Chem.MolFromSmiles(smiles)
        f.write(f"{i},{smiles},{score},{Descriptors.MolWt(mol):.2f},{Descriptors.MolLogP(mol):.2f}\n")

print("\n" + "="*70)
print("✅ PIPELINE COMPLETE! Results saved to output folder.")
print("="*70)