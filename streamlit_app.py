import streamlit as st
import requests
import py3Dmol

st.set_page_config(layout="wide")
st.title("🧪 3D Molecule Visualizer")

# Input
input_type = st.radio("Input Type", ("Molecule Name", "SMILES"))
user_input = st.text_input(f"Enter the {input_type}:")

# Convert molecule name to SMILES using PubChem
def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/TXT"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    return None

# Visualize molecule using py3Dmol
def show_3d_molecule(smiles):
    mol = py3Dmol.view(width=700, height=500)
    mol.addModel(smiles, "smi")
    mol.setStyle({'stick': {}})
    mol.setBackgroundColor("white")
    mol.zoomTo()
    return mol

# Display molecule
if user_input:
    if input_type == "Molecule Name":
        smiles = get_smiles_from_name(user_input)
        if not smiles:
            st.error("SMILES not found for the given molecule name.")
        else:
            st.success(f"Found SMILES: {smiles}")
    else:
        smiles = user_input

    if smiles:
        st.subheader("🧬 3D Structure")
        mol = show_3d_molecule(smiles)
        mol_html = mol._make_html()
        st.components.v1.html(mol_html, height=500, width=700)
        st.info("Use your mouse to rotate, zoom, and explore the structure.")

