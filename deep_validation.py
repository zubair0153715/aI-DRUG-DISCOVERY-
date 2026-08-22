"""
DEEP VALIDATION - PubChem Novelty + PAINS/Brenk + Synthetic Accessibility
"""

import requests
import urllib.parse
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Top candidates (from your results)
candidates = [
    "c1ccncc1CC(C)OCCO",
    "FCOCCOCCOF",
    "c1ccncc1OCCOOCCOC(=O)NCl",
    "CC(C)c1ccncc1COCCOC(=O)NCOC",
    "ClCOCCOCOCOCCO",
]

# ============================================
# 1. PUBCHEM NOVELTY CHECK
# ============================================
def pubchem_novelty_check(smiles):
    """Check if molecule exists in PubChem"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{urllib.parse.quote(smiles)}/cids/JSON"
    
    try:
        response = requests.get(url, timeout=10)
        data = response.json()
        
        if "IdentifierList" in data:
            cids = data["IdentifierList"].get("CID", [])
            return {
                "exists": True,
                "cids": cids,
                "novel": False,
                "message": f"Found in PubChem (CID: {cids[0] if cids else 'N/A'})"
            }
        else:
            return {
                "exists": False,
                "cids": [],
                "novel": True,
                "message": "NOVEL - Not found in PubChem!"
            }
    except Exception as e:
        return {
            "exists": False,
            "cids": [],
            "novel": None,
            "message": f"Check failed: {str(e)[:50]}"
        }

# ============================================
# 2. PAINS ALERTS
# ============================================
def pains_check(smiles):
    """Check for PAINS (Pan Assay Interference Compounds)"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"alerts": [], "pass": None}
    
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    
    entry = catalog.GetFirstMatch(mol)
    
    if entry:
        return {
            "alerts": [entry.GetDescription()],
            "pass": False,
            "message": f"⚠️ PAINS Alert: {entry.GetDescription()}"
        }
    else:
        return {
            "alerts": [],
            "pass": True,
            "message": "✅ PAINS: No alerts"
        }

# ============================================
# 3. BRENK ALERTS
# ============================================
def brenk_check(smiles):
    """Check for Brenk alerts (toxic/unstable groups)"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"alerts": [], "pass": None}
    
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog(params)
    
    entry = catalog.GetFirstMatch(mol)
    
    if entry:
        return {
            "alerts": [entry.GetDescription()],
            "pass": False,
            "message": f"⚠️ Brenk Alert: {entry.GetDescription()}"
        }
    else:
        return {
            "alerts": [],
            "pass": True,
            "message": "✅ Brenk: No alerts"
        }

# ============================================
# 4. SYNTHETIC ACCESSIBILITY (SA Score)
# ============================================
def synthetic_accessibility(smiles):
    """Calculate Synthetic Accessibility Score (1=easy, 10=hard)"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"score": None, "verdict": "N/A"}
    
    # Simplified SA Score calculation
    mw = Descriptors.MolWt(mol)
    rings = Descriptors.RingCount(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    chiral_centers = len(Chem.FindMolChiralCenters(mol))
    
    # Approximate SA score
    sa_score = 1 + (mw - 150) / 100 + rings * 0.3 + rot_bonds * 0.1 + chiral_centers * 0.5
    sa_score = min(10, max(1, sa_score))
    
    if sa_score < 3:
        verdict = "🟢 EASY to synthesize"
    elif sa_score < 5:
        verdict = "🟡 MODERATE"
    else:
        verdict = "🔴 DIFFICULT"
    
    return {
        "score": round(sa_score, 2),
        "verdict": verdict
    }

# ============================================
# RUN DEEP VALIDATION
# ============================================
print("="*80)
print("🔬 DEEP VALIDATION - PubChem + PAINS + Brenk + SA Score")
print("="*80)

for i, smiles in enumerate(candidates, 1):
    print(f"\n{'='*80}")
    print(f"CANDIDATE #{i}: {smiles}")
    print(f"{'='*80}")
    
    # PubChem
    print("\n📚 PubChem Novelty Check...")
    pubchem = pubchem_novelty_check(smiles)
    print(f"   {pubchem['message']}")
    
    # PAINS
    print("\n⚠️ PAINS Alert Check...")
    pains = pains_check(smiles)
    print(f"   {pains['message']}")
    
    # Brenk
    print("\n⚠️ Brenk Alert Check...")
    brenk = brenk_check(smiles)
    print(f"   {brenk['message']}")
    
    # Synthetic Accessibility
    print("\n🔧 Synthetic Accessibility...")
    sa = synthetic_accessibility(smiles)
    print(f"   SA Score: {sa['score']} - {sa['verdict']}")

print("\n" + "="*80)
print("✅ DEEP VALIDATION COMPLETE!")
print("="*80)