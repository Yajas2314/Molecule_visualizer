import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import requests
import numpy as np
import open3d as o3d
import tempfile
import os

st.title("Molecule Visualizer and .ply Exporter")

# Function to get SMILES from molecule name
def get_smiles_from_name(name):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{name}/smiles"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    else:
        return None

# Function to visualize molecule
def show_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    mb = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=600, height=400)
    viewer.addModel(mb, 'mol')
    viewer.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    viewer.addLabels()
    viewer.zoomTo()
    return viewer

# Function to generate colored .ply file
def generate_ply(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    points = []
    colors = []

    atom_colors = {
        'H': [1.0, 1.0, 1.0],
        'C': [0.5, 0.5, 0.5],
        'O': [1.0, 0.0, 0.0],
        'N': [0.0, 0.0, 1.0],
        'S': [1.0, 1.0, 0.0],
    }

    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        points.append([pos.x, pos.y, pos.z])
        sym = atom.GetSymbol()
        color = atom_colors.get(sym, [0.0, 1.0, 0.0])  # default green
        colors.append(color)

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(np.array(points))
    pcd.colors = o3d.utility.Vector3dVector(np.array(colors))

    ply_path = os.path.join(tempfile.gettempdir(), "molecule.ply")
    o3d.io.write_point_cloud(ply_path, pcd)

    return ply_path

# --- Streamlit Interface ---
user_input = st.text_input("Enter SMILES or Molecule Name")

if user_input:
    if all(c.isalpha() or c.isspace() for c in user_input):  # likely name
        smiles = get_smiles_from_name(user_input)
        if not smiles:
            st.error("Could not find SMILES for the given molecule name.")
        else:
            st.success(f"Found SMILES: {smiles}")
    else:
        smiles = user_input.strip()

    if smiles:
        viewer = show_molecule(smiles)
        st.components.v1.html(viewer._make_html(), height=450)

        if st.button("Generate .ply File"):
            try:
                ply_file = generate_ply(smiles)
                with open(ply_file, "rb") as f:
                    st.download_button("Download .ply", f, file_name="molecule.ply")
            except Exception as e:
                st.error(f"Error generating .ply: {str(e)}")
