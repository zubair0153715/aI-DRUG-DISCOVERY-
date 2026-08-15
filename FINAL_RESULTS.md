# AI Drug Discovery - Complete Results

## 3-Round Pipeline Results

### Round 1: Initial Generation + Docking
- 20 molecules generated
- Best binding: -7.44 kcal/mol (Compound 2)
- Target: EGFR (Cancer Protein)

### Round 2: ADMET Analysis (5 molecules)
- All passed Lipinski, PAINS = 0
- hERG: All safe
- Compound 3: Best DILI (0.207)

### Round 3: Optimization (18 variants)
- 10 variants + 8 balanced variants
- SAR Trend: OH > F > Cl > Br

## Final Lead Candidates

| Rank | SMILES | Binding | LogS | DILI | QED |
|------|--------|---------|------|------|-----|
| 1 | CC(CNC(=O)c1cccnc1)COCCCC(F)(F)F | -6.23 | -2.195 | 0.207 | 0.751 |
| 2 | COCCCOCC(C)CNC(=O)c1cccnc1 | -5.67 | -0.075 | 0.333 | 0.688 |
| 3 | COCCOCC(C)CNC(=O)c1cccnc1 | -5.54 | +0.048 | 0.260 | 0.704 |

## Key SAR Discovery
- Oxygen-rich tails: Better solubility, lower toxicity
- Fluorinated tails: Better binding, poor solubility
- Best balance: C5 (moderate both)

## Tools Used
- Python, RDKit, PyTorch, REINVENT 4
- AutoDock Vina, Meeko, Open Babel
- Google Colab, ADMETlab 2.0
- Biopython, OpenMM