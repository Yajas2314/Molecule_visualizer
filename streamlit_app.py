# streamlit_app.py - Pure Python Molecule Visualizer with Atom Labels and Enhanced Lone Pair Electrons

import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
from rdkit.Chem import AllChem
import py3Dmol
import requests
from PIL import Image
from io import BytesIO
import math

# Atomic data dictionary
atomic_data = {
    "H": {"number": 1, "mass": 1.008, "config": "1s1", "lone_pairs": 0},
    "C": {"number": 6, "mass": 12.011, "config": "1s2 2s2 2p2", "lone_pairs": 0},
    "N": {"number": 7, "mass": 14.007, "config": "1s2 2s2 2p3", "lone_pairs": 1},
    "O": {"number": 8, "mass": 15.999, "config": "1s2 2s2 2p4", "lone_pairs": 2},
    "F": {"number": 9, "mass": 18.998, "config": "1s2 2s2 2p5", "lone_pairs": 3},
    "Cl": {"number": 17, "mass": 35.45, "config": "[Ne] 3s2 3p5", "lone_pairs": 3},
    "Na": {"number": 11, "mass": 22.99, "config": "[Ne] 3s1", "lone_pairs": 0},
    "S": {"number": 16, "mass": 32.06, "config": "[Ne] 3s2 3p4", "lone_pairs": 2},
}

# Element color dictionary
element_colors = {
    "H": "white", "C": "black", "O": "red", "N": "blue", "S": "yellow",
    "Cl": "green", "F": "lime", "Na": "orange"
}

def fetch_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        response = requests.get(url)
        data = response.json()
        return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    except:
        return None

def render_2d_molecule(mol):
    img = Draw.MolToImage(mol, size=(300, 300))
    return img

def render_3d_view_html(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    mb = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(mb, 'mol')

    # Add atom labels as part of sphere style with specific labels
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        label = atom.GetSymbol()
        viewer.setStyle({"model": -1, "serial": idx + 1}, {
            "sphere": {"scale": 0.3, "color": element_colors.get(label, "gray")},
            "label": {
                "text": label,
                "font": "Arial",
                "scale": 0.4,
                "fontColor": "black",
                "showBackground": False
            }
        })

    # Add lone pair electrons
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        lp_count = atomic_data.get(symbol, {}).get("lone_pairs", 0)

        radius = 0.3  # Bigger size
        offset = 1.0  # Farther from nucleus

        for i in range(lp_count):
            angle = (i + 1) * (360 / (lp_count + 1))
            dx = offset * math.cos(math.radians(angle))
            dy = offset * math.sin(math.radians(angle))
            viewer.addSphere({
                'center': {'x': pos.x + dx, 'y': pos.y + dy, 'z': pos.z + 0.5},
                'radius': radius,
                'color': 'gray',
                'opacity': 0.5
            })

    viewer.zoomTo()
    return viewer

def display_element_info(mol):
    elements = set([atom.GetSymbol() for atom in mol.GetAtoms()])
    st.subheader("🧠 Element Info")
    for elem in sorted(elements):
        data = atomic_data.get(elem)
        if data:
            st.markdown(f"**{elem}**  →  Atomic No: {data['number']}, Mass: {data['mass']}, Lone Pairs: {data['lone_pairs']}, Configuration: `{data['config']}`")
        else:
            st.markdown(f"**{elem}**  →  No data available.")

st.set_page_config(page_title="Molecule Visualizer", layout="wide")
st.title("🧪 Molecule Visualizer with Atom Labels & Lone Pairs")
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
            st.image(render_2d_molecule(mol))

            st.subheader("📘 Molecular Info")
            st.write(f"**SMILES:** {smiles}")
            st.write(f"**Formula:** {rdMolDescriptors.CalcMolFormula(mol)}")
            st.write(f"**Molecular Weight:** {Descriptors.ExactMolWt(mol):.2f} g/mol")

            display_element_info(mol)

        with col2:
            st.subheader("🧪 3D Structure (With Labels & Lone Pairs)")
            viewer = render_3d_view_html(smiles)
            components.html(viewer._make_html(), height=400)
