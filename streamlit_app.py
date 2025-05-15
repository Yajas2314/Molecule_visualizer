import streamlit as st
import py3Dmol
import requests

st.set_page_config(page_title="3D Molecule Visualizer", layout="centered")

st.title("🧪 3D Molecule Visualizer with AR/VR")

# Input from user
input_option = st.radio("Select Input Type:", ("Molecule Name", "SMILES Code"))

if input_option == "Molecule Name":
    molecule_name = st.text_input("Enter Molecule Name (e.g., ethanol, benzene):")
else:
    smiles = st.text_input("Enter SMILES Code (e.g., CCO for ethanol):")

# Function to fetch SMILES from name using PubChem
def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
    try:
        res = requests.get(url)
        res.raise_for_status()
        data = res.json()
        return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    except Exception:
        return None

# Function to show 3D molecule using py3Dmol
def show_molecule(smiles_str):
    mol_view = py3Dmol.view(width=500, height=400)
    mol_view.addModel(smiles_str, "smi")
    mol_view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    mol_view.zoomTo()
    return mol_view

# Run if input is given
if (input_option == "Molecule Name" and molecule_name) or (input_option == "SMILES Code" and smiles):
    if input_option == "Molecule Name":
        smiles = get_smiles_from_name(molecule_name)
        if smiles:
            st.success(f"Found SMILES: {smiles}")
        else:
            st.error("Molecule not found. Please check the spelling.")
            st.stop()

    st.subheader("🔬 3D Structure")
    viewer = show_molecule(smiles)
    viewer_html = viewer._make_html()
    st.components.v1.html(viewer_html, height=400)
