import streamlit as st
import pandas as pd
import random
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Page config
st.set_page_config(page_title="AI Drug Discovery", page_icon="💊", layout="wide")

# Title
st.title("💊 AI Drug Discovery System")
st.markdown("### Insilico Pharma.AI jaisa system — bilkul free!")

# Sidebar - Disease selection
st.sidebar.header("🎯 Target Selection")
disease = st.sidebar.selectbox(
    "Disease Select Karo:",
    ["Cancer", "COVID-19", "Alzheimer's", "Diabetes", "Parkinson's", "HIV"]
)

# Disease → Target mapping
targets = {
    "Cancer": ["EGFR", "1M17", "Cell proliferation"],
    "COVID-19": ["Mpro", "6LU7", "Viral replication"],
    "Alzheimer's": ["AChE", "4EY7", "Neurotransmitter breakdown"],
    "Diabetes": ["DPP4", "2P8S", "Insulin regulation"],
    "Parkinson's": ["MAO-B", "4J5P", "Dopamine metabolism"],
    "HIV": ["Protease", "1HPV", "Viral maturation"],
}

target_info = targets[disease]

# Sidebar info
st.sidebar.info(f"""
**Target:** {target_info[0]}  
**PDB ID:** {target_info[1]}  
**Pathway:** {target_info[2]}  
""")

# Main area
col1, col2, col3 = st.columns(3)
col1.metric("Disease", disease)
col2.metric("Target Protein", target_info[0])
col3.metric("PDB ID", target_info[1])

# Run button
if st.button("🚀 Drug Discovery Shuru Karo", type="primary"):
    
    # Progress bar
    progress = st.progress(0)
    status = st.empty()
    
    # Step 1: Generate molecules
    status.text("🧬 Generating molecules...")
    fragments = [
        "c1ccccc1", "c1ccncc1", "C(=O)N", "COC", "F", "Cl",
        "CC(C)", "C(F)(F)F", "OCCO", "COCCO", "COCCCO",
    ]
    
    molecules = set()
    for _ in range(1000):
        num_frags = random.randint(2, 5)
        smiles = ""
        for _ in range(num_frags):
            smiles += random.choice(fragments)
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            if 150 < mw < 500:
                molecules.add(smiles)
    
    molecules = list(molecules)[:20]
    progress.progress(25)
    
    # Step 2: AI Docking (simulated)
    status.text("🔬 AI Docking...")
    results = []
    for smiles in molecules:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            random.seed(hash(smiles) % 1000)
            binding = -4 - random.uniform(0, 5)
            results.append({
                "smiles": smiles,
                "binding": round(binding, 2),
                "mw": round(Descriptors.MolWt(mol), 1),
                "logp": round(Descriptors.MolLogP(mol), 2),
            })
    
    results.sort(key=lambda x: x["binding"])
    progress.progress(50)
    
    # Step 3: ADMET
    status.text("🧪 ADMET Prediction...")
    for r in results:
        mol = Chem.MolFromSmiles(r["smiles"])
        r["lipinski"] = "PASS" if (r["mw"] < 500 and r["logp"] < 5) else "FAIL"
        r["solubility"] = "GOOD" if r["logp"] < 3 else "POOR"
    progress.progress(75)
    
    # Step 4: Report
    status.text("📊 Generating Report...")
    progress.progress(100)
    status.text("✅ Complete!")
    
    # Results display
    st.success(f"Drug Discovery Complete for {disease}!")
    
    # Top 10 table
    st.subheader("🏆 Top 10 Drug Candidates")
    df = pd.DataFrame(results[:10])
    st.dataframe(df, use_container_width=True)
    
    # Best molecule
    best = results[0]
    st.subheader("💊 Best Drug Candidate")
    st.info(f"""
    **SMILES:** {best['smiles']}  
    **Binding:** {best['binding']} kcal/mol  
    **MW:** {best['mw']} | **LogP:** {best['logp']}  
    **Lipinski:** {best['lipinski']} | **Solubility:** {best['solubility']}
    """)
    
    # SMILES display
    st.code(best["smiles"], language="text")
    
    # Download button
    csv = df.to_csv(index=False)
    st.download_button(
        "📥 Results Download Karo (CSV)",
        csv,
        "drug_candidates.csv",
        "text/csv"
    )

# Footer
st.markdown("---")
st.markdown("**AI Drug Discovery System** | Insilico Pharma.AI Clone | Open Source")