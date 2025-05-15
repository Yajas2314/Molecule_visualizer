import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import requests

# Title
st.set_page_config(layout="wide")
st.markdown("<h1 style='text-align: center;'>🧪 3D Molecule Visualizer</h1>", unsafe_allow_html=True)
st.write("Enter either a **molecule name** or a **SMILES code** to view its 3D structure.")

# Sidebar Input
st.sidebar.header("Input Molecule")
mol_name = st.sidebar.text_input("Enter Molecule Name (e.g. glucose):")
smiles_input = st.sidebar.text_input("Or Enter SMILES Code (e.g. C(CO)O):")

style = st.sidebar.selectbox("Choose Visualization Style:", ['stick', 'sphere', 'line', 'cartoon'])

# Convert molecule name to SMILES
def get_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        response = requests.get(url)
        data = response.json()
        return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    except:
        return None

# Generate 3D structure from SMILES
def show_molecule(smiles, style="stick"):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    mol_block = Chem.MolToMolBlock(mol)

    viewer = py3Dmol.view(width=600, height=400)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({style: {}})
    viewer.zoomTo()
    return viewer

# Main Logic
if smiles_input:
    smiles = smiles_input
elif mol_name:
    smiles = get_smiles_from_name(mol_name)
else:
    smiles = None

if smiles:
    st.markdown("### 🔬 3D Structure Below:")
    viewer = show_molecule(smiles, style)
    viewer_html = viewer._make_html()
    st.components.v1.html(viewer_html, height=500, width=700)
else:
    st.warning("Please enter a molecule name or SMILES code.")
