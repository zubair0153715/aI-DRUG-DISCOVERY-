# 🧬 AI Drug Discovery System

A complete AI-powered drug discovery pipeline built with open-source tools. This system generates novel drug-like molecules, validates them, and performs molecular docking against cancer protein EGFR.

## 🎯 Target
- **Protein:** EGFR (Epidermal Growth Factor Receptor)
- **PDB ID:** 1M17
- **Disease:** Cancer
- **Binding Pocket Center:** (21.96, 6.05, 57.02)

## 🏆 Best Drug Candidate
SMILES: CN(C)c1ccccc1C(=O)NS(=O)(=O)Nc1ccccc1F
Binding Affinity: -7.44 kcal/mol (AutoDock Vina)
Molecular Weight: 337.38 g/mol
LogP: 1.98
TPSA: 78.51

## 📊 Docking Results (AutoDock Vina)

| Rank | SMILES | Binding (kcal/mol) | Verdict |
|------|--------|-------------------|---------|
| 1 | CN(C)c1ccccc1C(=O)NS(=O)(=O)Nc1ccccc1F | -7.44 | GOOD BINDER |
| 2 | OS(=O)(=O)NNCCc1ccncc1c1ccccc1Cl | -6.84 | GOOD BINDER |
| 3 | c1ccncc1C(=O)NCC(C)COCCCC(F)(F)F | -6.46 | GOOD BINDER |
| 4 | COCc1ccccc1OS(=O)(=O)NCC(C)C(=O)NNCC | -6.12 | GOOD BINDER |
| 5 | OCc1ccccc1Cl | -4.74 | WEAK BINDER |

## 🔬 ADMET Analysis Results

All 5 molecules passed ADMET screening:
- ✅ Molecular Weight: All < 500 g/mol
- ✅ LogP: All within 0-5 range
- ✅ TPSA: All < 140 (good absorption)
- ✅ H-Bond Donors: All ≤ 5
- ✅ H-Bond Acceptors: All ≤ 10
- ✅ No major toxicity flags

## 🔧 System Components

| Component | Technology |
|-----------|-----------|
| AI Model | REINVENT 4 (AstraZeneca) |
| Chemistry Engine | RDKit |
| Deep Learning | PyTorch |
| Docking | AutoDock Vina |
| Protein Analysis | Biopython + OpenMM |
| Ligand Preparation | Meeko |
| Cloud Computing | Google Colab |

## 📁 Output Files

- `generated_molecules.csv` — 20 AI-generated molecules
- `top_drug_candidates.csv` — Top 5 candidates
- `final_docking_results.csv` — Real docking scores
- `admet_results.csv` — ADMET analysis
- `ai_drugs.png` — Molecule visualizations

## 🚀 How to Run

### Local Setup
python -m venv drug_ai
drug_ai\Scripts\activate
pip install rdkit numpy pandas scikit-learn matplotlib biopython torch meeko
python train_simple_model.py
python generate_from_model.py
python admet_analysis.py

### Docking (Google Colab)
1. Upload docking scripts to Google Colab
2. Install AutoDock Vina
3. Run docking against EGFR protein
4. Get binding affinities

## 📈 Pipeline Flow

SMILES Input -> AI Model Training -> Molecule Generation -> Drug-likeness Filter -> Protein Download -> Ligand Preparation -> Molecular Docking -> Binding Analysis -> ADMET Screening -> Final Drug Candidates

## 💡 Key Features

- Zero-cost setup (all open-source tools)
- Cloud-based docking for real results
- Automated molecule generation
- Drug-likeness validation (Lipinski's Rule)
- ADMET safety screening
- Professional CSV outputs

## 📧 Contact

Open to collaboration, research, and freelance opportunities.

## 📄 License

MIT License