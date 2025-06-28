# streamlit_app.py - Pure Python Molecule Visualizer

# streamlit_app.py - Pure Python Molecule Visualizer with HTML-based 3D Viewer

import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
from rdkit.Chem import AllChem
import py3Dmol
import requests
from PIL import Image
from io import BytesIO

# Atomic data dictionary
atomic_data = {
    "H": {"number": 1, "mass": 1.008, "config": "1s1"},
    "C": {"number": 6, "mass": 12.011, "config": "1s2 2s2 2p2"},
    "N": {"number": 7, "mass": 14.007, "config": "1s2 2s2 2p3"},
    "O": {"number": 8, "mass": 15.999, "config": "1s2 2s2 2p4"},
    "F": {"number": 9, "mass": 18.998, "config": "1s2 2s2 2p5"},
    "Cl": {"number": 17, "mass": 35.45, "config": "[Ne] 3s2 3p5"},
    "Na": {"number": 11, "mass": 22.99, "config": "[Ne] 3s1"},
    "S": {"number": 16, "mass": 32.06, "config": "[Ne] 3s2 3p4"},
    # Extend this dictionary as needed
}

# Element color dictionary for 3Dmol.js
element_colors = {
    "H": "white", "C": "black", "O": "red", "N": "blue", "S": "yellow",
    "Cl": "green", "F": "lime", "Na": "orange"
}

# Get SMILES from PubChem if user input is a name
def fetch_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        response = requests.get(url)
        data = response.json()
        return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    except:
        return None

# Render 2D image from RDKit molecule
def render_2d_molecule(mol):
    img = Draw.MolToImage(mol, size=(300, 300))
    return img

# Render 3D viewer using py3Dmol and display via HTML iframe

def render_3d_view_html(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    mb = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(mb, 'mol')
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
    viewer.zoomTo()
    return viewer

# Show element-wise data below structure
def display_element_info(mol):
    elements = set([atom.GetSymbol() for atom in mol.GetAtoms()])
    st.subheader("🧠 Element Info")
    for elem in sorted(elements):
        data = atomic_data.get(elem)
        if data:
            st.markdown(f"**{elem}**  →  Atomic No: {data['number']}, Mass: {data['mass']}, Configuration: `{data['config']}`")
        else:
            st.markdown(f"**{elem}**  →  No data available.")

# App UI starts here
st.set_page_config(page_title="Molecule Visualizer", layout="wide")
st.title("🧪 Molecule Visualizer with Atom Data")
user_input = st.text_input("Enter molecule name or SMILES (e.g. Ethanol or CCO):")

if user_input:
    smiles = user_input if Chem.MolFromSmiles(user_input) else fetch_smiles_from_name(user_input)

    if not smiles:
        st.error("Could not resolve input to a valid SMILES string.")
    else:
        mol = Chem.MolFromSmiles(smiles)
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("🧬 2D Structure")
            img = render_2d_molecule(mol)
            st.image(img)

            st.subheader("📘 Molecular Info")
            st.write(f"**SMILES:** {smiles}")
            st.write(f"**Formula:** {rdMolDescriptors.CalcMolFormula(mol)}")
            st.write(f"**Molecular Weight:** {Descriptors.ExactMolWt(mol):.2f} g/mol")

            display_element_info(mol)

        with col2:
            st.subheader("🧪 3D Structure (Color by Element)")
            viewer = render_3d_view_html(smiles)
            components.html(viewer._make_html(), height=400)
