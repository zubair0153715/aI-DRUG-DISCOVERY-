import os
from rdkit import Chem
from rdkit.Chem import Descriptors

print("="*70)
print("🔬 ADMET ANALYSIS - AI GENERATED DRUG CANDIDATES")
print("="*70)

# Aapke 5 best molecules
molecules = [
    ("Molecule 1 (Best)", "CN(C)c1ccccc1C(=O)NS(=O)(=O)Nc1ccccc1F"),
    ("Molecule 2", "OS(=O)(=O)NNCCc1ccncc1c1ccccc1Cl"),
    ("Molecule 3", "c1ccncc1C(=O)NCC(C)COCCCC(F)(F)F"),
    ("Molecule 4", "COCc1ccccc1OS(=O)(=O)NCC(C)C(=O)NNCC"),
    ("Molecule 5", "OCc1ccccc1Cl"),
]

print("\n📊 KEY PROPERTIES FOR DRUG SAFETY:")
print("-"*70)

for name, smiles in molecules:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rings = Descriptors.RingCount(mol)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        
        print(f"\n🧪 {name}")
        print(f"   SMILES: {smiles}")
        print(f"   Molecular Weight: {mw:.2f} g/mol (Ideal: <500)")
        print(f"   LogP: {logp:.2f} (Ideal: 0-5)")
        print(f"   H-Bond Donors: {hbd} (Ideal: ≤5)")
        print(f"   H-Bond Acceptors: {hba} (Ideal: ≤10)")
        print(f"   TPSA: {tpsa:.2f} (Ideal: <140)")
        print(f"   Rings: {rings}")
        print(f"   Rotatable Bonds: {rot_bonds}")
        
        # Safety assessment
        issues = []
        if mw > 500: issues.append("High MW")
        if logp > 5: issues.append("Too lipophilic")
        if hbd > 5: issues.append("Too many HBD")
        if hba > 10: issues.append("Too many HBA")
        if tpsa > 140: issues.append("High TPSA (poor absorption)")
        
        if issues:
            print(f"   ⚠️ Issues: {', '.join(issues)}")
        else:
            print(f"   ✅ No major ADMET issues detected")

# Save results
output_file = os.path.join(os.getcwd(), "output", "admet_results.csv")
with open(output_file, "w") as f:
    f.write("Molecule,SMILES,MW,LogP,HBD,HBA,TPSA,Rings,RotBonds\n")
    for name, smiles in molecules:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            f.write(f"{name},{smiles},{Descriptors.MolWt(mol):.2f},{Descriptors.MolLogP(mol):.2f},{Descriptors.NumHDonors(mol)},{Descriptors.NumHAcceptors(mol)},{Descriptors.TPSA(mol):.2f},{Descriptors.RingCount(mol)},{Descriptors.NumRotatableBonds(mol)}\n")

print("\n" + "="*70)
print("✅ ADMET Analysis Complete! Results saved to output folder.")
print("="*70)