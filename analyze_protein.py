import os
from Bio.PDB import PDBParser

print("="*60)
print("PROTEIN TARGET ANALYSIS...")
print("="*60)

# Protein file ka path
pdb_file = os.path.join(os.getcwd(), "data", "1m17.pdb")

if os.path.exists(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", pdb_file)
    
    print(f"✅ Protein Loaded: {pdb_file}")
    print(f"\n📊 PROTEIN INFORMATION:")
    print("-"*60)
    
    for model in structure:
        for chain in model:
            residues = list(chain.get_residues())
            print(f"Chain ID: {chain.id}")
            print(f"Total Residues (Amino Acids): {len(residues)}")
            print(f"First 5 Residues: {[res.resname for res in residues[:5]]}")
            break
        break
    
    print("\n✅ Target Protein Ready for Docking!")
else:
    print("Protein file not found. Trying COVID protein...")
    pdb_file = os.path.join(os.getcwd(), "data", "6lu7.pdb")
    if os.path.exists(pdb_file):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("target", pdb_file)
        print(f"✅ COVID Protein Loaded: {pdb_file}")
        print("✅ Target Ready for Docking!")
    else:
        print("❌ No protein found. Check download.")

print("="*60)