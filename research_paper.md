Title
Open-Source AI-Driven Drug Discovery Pipeline: Generation, Docking, and ADMET Optimization of Novel EGFR Inhibitors

Authors
Muhammad Zubair
Independent Researcher
Email: zubair0153715@gmail.com
GitHub: https://github.com/zubair0153715/aI-DRUG-DISCOVERY-

Abstract
Drug discovery is a time-intensive and costly process, often requiring over a decade and billions of dollars to bring a single drug to market. Recent advances in artificial intelligence (AI) and cloud computing have opened new possibilities for accelerating early-stage drug discovery. In this study, we present a complete open-source AI-driven drug discovery pipeline that integrates molecular generation, molecular docking, and ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) prediction. Using REINVENT 4, RDKit, PyTorch, AutoDock Vina, and ADMETlab 2.0, we generated 20 novel drug-like molecules targeting the Epidermal Growth Factor Receptor (EGFR, PDB ID: 1M17), a validated cancer target. Initial docking identified a lead compound with a binding affinity of -7.44 kcal/mol. Subsequent optimization through 18 structural variants revealed a clear structure-activity relationship (SAR), demonstrating that oxygen-rich tails significantly improve solubility and reduce toxicity, while fluorinated tails enhance binding affinity. Our best-balanced lead candidates achieved binding affinities between -5.54 and -6.23 kcal/mol with favorable ADMET profiles. This work demonstrates that open-source tools can enable accessible, low-cost computational drug discovery, providing a foundation for further experimental validation.

Keywords: AI drug discovery, EGFR, molecular docking, ADMET, REINVENT 4, AutoDock Vina, RDKit, open-source

1. Introduction
The traditional drug discovery process is characterized by high failure rates, with approximately 90% of candidates failing clinical trials. The average cost to develop a new drug exceeds $2.6 billion, and the timeline spans 10–15 years. These challenges necessitate innovative approaches to accelerate early-stage discovery.

Artificial intelligence has emerged as a transformative force in drug discovery. Machine learning models can now generate novel molecules, predict their binding to target proteins, and assess their safety profiles—all before entering a laboratory. Notable successes include Insilico Medicine's AI-generated drug candidate for idiopathic pulmonary fibrosis, which entered clinical trials in under 18 months.

This study explores the feasibility of building a complete AI drug discovery pipeline using only open-source tools and free cloud resources. We target EGFR, a well-established receptor tyrosine kinase implicated in non-small cell lung cancer (NSCLC), colorectal cancer, and glioblastoma. EGFR inhibitors such as gefitinib and erlotinib have shown clinical efficacy, making EGFR an ideal target for demonstrating our pipeline.

2. Methods
2.1 System Architecture
Our pipeline consists of four major components:

Molecular Generation: REINVENT 4 (AstraZeneca) and custom RDKit-based models generate drug-like molecules.

Drug-likeness Filtering: Lipinski's Rule of Five and QED scores filter compounds.

Molecular Docking: AutoDock Vina predicts binding affinities against EGFR.

ADMET Prediction: ADMETlab 2.0 evaluates 89 pharmacokinetic and toxicity parameters.

2.2 Molecular Generation
We utilized REINVENT 4, an open-source framework from MolecularAI (AstraZeneca), for molecular generation. A custom prior model was trained using 20 chemical fragments derived from common drug structures. The model generated 20 novel molecules with drug-like properties.

2.3 Target Protein Preparation
The EGFR crystal structure (PDB ID: 1M17, 333 amino acid residues) was downloaded from the RCSB Protein Data Bank. The protein was prepared by removing water molecules and heteroatoms, retaining only ATOM records for docking. The binding pocket center was determined from HETATM records of the co-crystallized ligand, yielding coordinates (X: 21.96, Y: 6.05, Z: 57.02).

2.4 Molecular Docking
AutoDock Vina was used for molecular docking on Google Colab (free cloud platform). Ligands were prepared using RDKit and Open Babel for conversion to PDBQT format. Docking was performed with exhaustiveness of 8 and a grid box of 20 × 20 × 20 Å centered at the binding pocket.

2.5 ADMET Analysis
ADMET predictions were performed using ADMETlab 2.0. Each compound was evaluated for 89 parameters, including:

Physicochemical properties (MW, LogP, LogS, TPSA)

Absorption (HIA, Caco-2, F20)

Distribution (BBB, PPB)

Metabolism (CYP450 inhibition)

Toxicity (DILI, hERG, Ames, Carcinogenicity)

2.6 Optimization Strategy
Three rounds of optimization were performed:

Round 1: Initial 5 molecules analyzed

Round 2: 10 structural variants generated

Round 3: 8 balanced variants combining binding and safety features

3. Results
3.1 Molecular Generation
The AI model successfully generated 20 novel drug-like molecules. All 20 molecules passed Lipinski's Rule of Five and showed PAINS = 0 (no structural alerts). Molecular weights ranged from 101 to 345 Da, and LogP values ranged from -1.79 to 4.32.

3.2 Initial Docking Results
From the initial set, 5 molecules were selected for docking against EGFR. The best binding affinity was -7.44 kcal/mol for Compound 2 (CN(C)c1ccccc1C(=O)NS(=O)(=O)Nc1ccccc1F).

Table 1: Initial Docking Results

Rank	SMILES	Binding (kcal/mol)
1	CN(C)c1ccccc1C(=O)NS(=O)(=O)Nc1ccccc1F	-7.44
2	OS(=O)(=O)NNCCc1ccncc1c1ccccc1Cl	-6.84
3	c1ccncc1C(=O)NCC(C)COCCCC(F)(F)F	-6.46
4	COCc1ccccc1OS(=O)(=O)NCC(C)C(=O)NNCC	-6.12
5	OCc1ccccc1Cl	-4.74
3.3 ADMET Analysis (Round 1)
Table 2: ADMET Round 1 Results

Compound	DILI	hERG	QED	Lipinski	Verdict
1	0.997	0.003	0.874	Pass	Reject (DILI)
2	0.947	0.090	0.428	Pass	Reject (DILI)
3	0.207	0.107	0.751	Pass	Modify
4	0.993	0.007	0.528	Pass	Reject (DILI)
5	0.434	0.026	0.635	Pass	Modify
3.4 Optimization Results
Table 3: Key ADMET Improvements After Optimization

Property	Original (C1)	Best Optimized (C4)	Improvement
LogS	-2.195	+0.048	100x better
PPB	68.98%	13.38%	5x better
Free Fraction	28.77%	75.06%	2.6x better
CYP3A4	0.377	0.048	8x better
BCF	0.862	0.150	6x better
DILI	0.207	0.260	Similar
3.5 Final Lead Candidates
Table 4: Final Lead Candidates

Rank	SMILES	Binding (kcal/mol)	LogS	DILI	QED
1	CC(CNC(=O)c1cccnc1)COCCCC(F)(F)F	-6.23	-2.195	0.207	0.751
2	COCCCOCC(C)CNC(=O)c1cccnc1	-5.67	-0.075	0.333	0.688
3	COCCOCC(C)CNC(=O)c1cccnc1	-5.54	+0.048	0.260	0.704
4. Discussion
4.1 Structure-Activity Relationship (SAR)
A clear SAR trend emerged from our optimization study:

Fluorinated tails (CF3, CF2, CF): Higher binding affinity but poor solubility (LogS -1.5 to -2.8) and higher toxicity

Oxygen-rich tails (OH, ether): Excellent solubility (LogS +0.048) and lower toxicity, but reduced binding

This trade-off represents a classic medicinal chemistry challenge: balancing potency with developability.

4.2 Comparison with Existing Drugs
Approved EGFR inhibitors (gefitinib, erlotinib) have binding affinities ranging from -8.0 to -12.0 kcal/mol. Our best lead at -6.23 kcal/mol is a promising starting point requiring further optimization to reach clinical-grade potency.

4.3 Advantages of Open-Source Approach
Our pipeline demonstrates several advantages:

Cost: Complete pipeline runs on free tools

Accessibility: No specialized hardware required

Reproducibility: All code available on GitHub

Scalability: Can target any protein in PDB

4.4 Limitations
Docking scores are predictions, not experimental validations

ADMET predictions use machine learning models with inherent uncertainties

Permeability remains a challenge despite good solubility

5. Conclusion
We successfully built and validated a complete AI-driven drug discovery pipeline using open-source tools. The pipeline generated 20 novel drug-like molecules, identified lead compounds through molecular docking, and optimized them through ADMET screening. The best-balanced lead candidate (COCCOCC(C)CNC(=O)c1cccnc1) showed excellent solubility (LogS +0.048), high free fraction (75%), and low toxicity, though with moderate binding affinity (-5.54 kcal/mol). This work demonstrates that accessible computational drug discovery is feasible and provides a foundation for experimental validation.

6. Future Work
Molecular Dynamics simulations (100 ns) for lead compounds

Permeability optimization while maintaining solubility

Testing against additional cancer targets (HER2, BRAF)

Wet lab validation through collaboration with academic institutions

Integration of generative AI for targeted optimization

7. References
Trott, O., & Olson, A. J. (2010). AutoDock Vina: Improving the speed and accuracy of docking. Journal of Computational Chemistry, 31(2), 455-461.

Eberhardt, J., et al. (2021). AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. Journal of Chemical Information and Modeling.

Blaschke, T., et al. (2020). REINVENT 2.0: An AI Tool for De Novo Drug Design. Journal of Chemical Information and Modeling.

Lipinski, C. A. (2004). Lead- and drug-like compounds: The rule-of-five revolution. Drug Discovery Today: Technologies.

Xiong, G., et al. (2021). ADMETlab 2.0: An integrated online platform for accurate and comprehensive predictions of ADMET properties. Nucleic Acids Research.

Berman, H. M., et al. (2000). The Protein Data Bank. Nucleic Acids Research, 28(1), 235-242.

RDKit: Open-source cheminformatics. https://www.rdkit.org

MolecularAI. REINVENT 4. https://github.com/MolecularAI/REINVENT4

8. Supplementary Information
All code, data, and results are available at:
https://github.com/zubair0153715/aI-DRUG-DISCOVERY-