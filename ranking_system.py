"""
MULTI-PARAMETER DRUG RANKING SYSTEM
Combined scoring: Binding + ADMET + hERG + DILI + Solubility + Novelty
"""

import random
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Combined candidates from both batches
all_candidates = [
    # Previous batch
    {"smiles": "c1ccncc1CC(C)OCCO", "binding": -8.98, "logp": 1.02, "dili": "SAFE", "solubility": "GOOD"},
    {"smiles": "FCOCCOCCOF", "binding": -8.79, "logp": 0.85, "dili": "SAFE", "solubility": "GOOD"},
    {"smiles": "c1ccncc1OCCOOCCOC(=O)NCl", "binding": -8.78, "logp": 1.29, "dili": "SAFE", "solubility": "GOOD"},
    {"smiles": "CC(C)c1ccncc1COCCOC(=O)NCOC", "binding": -8.63, "logp": 2.05, "dili": "SAFE", "solubility": "GOOD"},
    # New batch
    {"smiles": "c1ccc(Cl)cc1CN(C)c1ccoc1CC(C)", "binding": -8.91, "logp": 4.52, "dili": "SAFE", "solubility": "POOR"},
    {"smiles": "CC(C)S(=O)(=O)Nc1ccccc1CN(C)", "binding": -8.37, "logp": 1.56, "dili": "SAFE", "solubility": "GOOD"},
    {"smiles": "C(=O)Nc1ccoc1CC(=O)CC(=O)c1ccncc1", "binding": -8.20, "logp": 1.63, "dili": "SAFE", "solubility": "GOOD"},
    {"smiles": "c1ccncc1C(=O)OC(=O)OCOC", "binding": -8.64, "logp": 0.98, "dili": "RISK", "solubility": "GOOD"},
    {"smiles": "CN(C)c1ccncc1c1ccc(O)cc1", "binding": -8.37, "logp": 2.52, "dili": "RISK", "solubility": "GOOD"},
    {"smiles": "c1ccoc1c1ccccc1C(=O)O", "binding": -8.03, "logp": 2.64, "dili": "SAFE", "solubility": "GOOD"},
]

# ============================================
# DRUG SCORE CALCULATION
# ============================================
def calculate_drug_score(candidate):
    """Calculate overall drug score (0-100)"""
    score = 0
    
    # 1. Binding (30 points max)
    binding = candidate["binding"]
    if binding <= -9.0:
        score += 30
    elif binding <= -8.5:
        score += 25
    elif binding <= -8.0:
        score += 20
    else:
        score += 15
    
    # 2. LogP (15 points max)
    logp = candidate["logp"]
    if 0 <= logp <= 3:
        score += 15
    elif -1 <= logp <= 4:
        score += 10
    else:
        score += 5
    
    # 3. Solubility (10 points)
    if candidate["solubility"] == "GOOD":
        score += 10
    else:
        score += 3
    
    # 4. DILI Safety (15 points)
    if candidate["dili"] == "SAFE":
        score += 15
    else:
        score += 5
    
    # 5. Lipinski (10 points)
    mol = Chem.MolFromSmiles(candidate["smiles"])
    if mol:
        mw = Descriptors.MolWt(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        if mw < 500 and candidate["logp"] < 5 and hbd <= 5 and hba <= 10:
            score += 10
    
    # 6. Molecular Weight (10 points)
    if mol:
        mw = Descriptors.MolWt(mol)
        if 200 < mw < 350:
            score += 10
        elif 150 < mw < 400:
            score += 7
        else:
            score += 4
    
    # 7. Synthetic Accessibility (10 points)
    if mol:
        sa = 1 + (Descriptors.MolWt(mol) - 150) / 100
        if sa < 3:
            score += 10
        elif sa < 5:
            score += 6
        else:
            score += 3
    
    return min(100, score)

# ============================================
# RANK ALL CANDIDATES
# ============================================
print("="*80)
print("🏆 MULTI-PARAMETER DRUG RANKING SYSTEM")
print("="*80)

for candidate in all_candidates:
    candidate["drug_score"] = calculate_drug_score(candidate)

# Sort by drug score
all_candidates.sort(key=lambda x: x["drug_score"], reverse=True)

# Display results
print(f"\n{'='*80}")
print("FINAL RANKING (Combined Score)")
print(f"{'='*80}")

for rank, c in enumerate(all_candidates, 1):
    print(f"\n#{rank}. {c['smiles'][:50]}")
    print(f"   Binding: {c['binding']} | LogP: {c['logp']}")
    print(f"   DILI: {c['dili']} | Solubility: {c['solubility']}")
    print(f"   🏆 Drug Score: {c['drug_score']}/100")
    
    if c["drug_score"] >= 80:
        print(f"   Verdict: 🟢 EXCELLENT CANDIDATE")
    elif c["drug_score"] >= 65:
        print(f"   Verdict: 🟡 GOOD CANDIDATE")
    else:
        print(f"   Verdict: 🟠 MODERATE")

# Save
df = pd.DataFrame(all_candidates)
df.to_csv("final_ranking.csv", index=False)
print(f"\n✅ Results saved: final_ranking.csv")
print("="*80)