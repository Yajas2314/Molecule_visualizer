import streamlit as st
import py3Dmol
from urllib.parse import quote
import requests

def get_smiles_from_name(name):
    """Fetch SMILES code from molecule name using PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(name)}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            data = response.json()
            smiles = data['PropertyTable']['Properties'][0]['IsomericSMILES']
            return smiles
        except Exception:
            return None
    return None

def show_molecule(smiles):
    """Render 3D molecule using py3Dmol."""
    mol = py3Dmol.view(width=400, height=400)
    mol.addModel(smiles, "smi")
    mol.setStyle({'stick': {}})
    mol.zoomTo()
    return mol

st.title("🧪 Molecule Visualizer")
input_type = st.radio("Choose input method:", ("Molecule Name", "SMILES Code"))

if input_type == "Molecule Name":
    name = st.text_input("Enter molecule name:")
    if name:
        smiles = get_smiles_from_name(name)
        if smiles:
            st.success(f"Found SMILES: `{smiles}`")
            st.subheader("3D Structure")
            mol = show_molecule(smiles)
            mol_html = mol._make_html()
            st.components.v1.html(mol_html, height=400)
        else:
            st.error("Could not fetch SMILES for the given molecule name.")
else:
    smiles = st.text_input("Enter SMILES code:")
    if smiles:
        st.subheader("3D Structure")
        mol = show_molecule(smiles)
        mol_html = mol._make_html()
        st.components.v1.html(mol_html, height=400)
