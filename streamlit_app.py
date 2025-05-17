import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import requests


def get_smiles_from_name(name):
    """Get SMILES from a molecule name using PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data['PropertyTable']['Properties'][0]['IsomericSMILES']
    return None


def draw_3d_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    mol_block = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=600, height=500)
    viewer.addModel(mol_block, 'mol')
    viewer.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    viewer.setBackgroundColor('white')
    viewer.zoomTo()
    return viewer._make_html()


st.set_page_config(page_title="Molecule Visualizer", layout="centered")
st.title("🧪 3D Molecule Visualizer with SMILES or Name")

user_input = st.text_input("Enter SMILES or Molecule Name:")

if user_input:
    if any(ch in user_input for ch in "=#()/\\"):  # Likely SMILES
        smiles = user_input
    else:
        smiles = get_smiles_from_name(user_input)

    if smiles:
        st.success(f"SMILES: `{smiles}`")
        html_view = draw_3d_molecule(smiles)
        st.components.v1.html(html_view, height=500, width=700)
    else:
        st.error("❌ Could not get SMILES for the given input.")
