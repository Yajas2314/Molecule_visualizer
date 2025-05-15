import streamlit as st
import py3Dmol
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

# Get SMILES from PubChem using molecule name
def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/TXT"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    return None

# Detect bond type (very basic, can be expanded)
def classify_compound(smiles):
    ionic_keywords = ['[Na+]', '[Cl-]', '[K+]', '[Ca+2]', '[Mg+2]']
    for ion in ionic_keywords:
        if ion in smiles:
            return "Likely Ionic"
    return "Likely Covalent"

# Render molecule in 3D using py3Dmol
def show_molecule(smiles, style='stick'):
    view = py3Dmol.view(width=500, height=400)
    view.addModel(smiles, 'smi')
    view.setStyle({style: {}})
    view.zoomTo()
    return view

# Streamlit app
def main():
    st.set_page_config(page_title="3D Molecule Visualizer", layout="wide")
    st.markdown("<h1 style='text-align: center;'>🧪 3D Molecule Visualizer</h1>", unsafe_allow_html=True)

    # Sidebar inputs
    with st.sidebar:
        st.header("Input Molecule")
        mol_name = st.text_input("Enter Molecule Name (e.g. glucose):")
        mol_smiles = st.text_input("Or Enter SMILES Code (e.g. C(CO)O):")

        style = st.selectbox("Choose Visualization Style:", ["stick", "line", "sphere", "cartoon", "cross"])

    # Determine which SMILES to use
    smiles = mol_smiles.strip()
    if not smiles and mol_name:
        smiles = get_smiles_from_name(mol_name.strip())

    if smiles:
        compound_type = classify_compound(smiles)

        st.subheader("🧬 3D Structure Viewer")
        col1, col2 = st.columns([1, 2])

        with col1:
            st.markdown(f"**Compound Type:** {compound_type}")
            st.markdown(f"**SMILES Code:** `{smiles}`")

        with col2:
            mol_view = show_molecule(smiles, style)
            mol_html = mol_view._make_html()
            st.components.v1.html(mol_html, height=420, width=500)

    else:
        st.warning("Please enter a valid molecule name or SMILES code.")

if __name__ == "__main__":
    main()
