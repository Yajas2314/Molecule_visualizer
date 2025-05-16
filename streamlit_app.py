import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import pyvista as pv
import numpy as np
import tempfile
import os

st.title("🔬 Molecule 3D Viewer and .PLY Exporter")
st.write("Enter a molecule name or SMILES code:")

# Convert molecule name to SMILES using CACTUS
def get_smiles_from_name(name):
    try:
        url = f"https://cactus.nci.nih.gov/chemical/structure/{name}/smiles"
        response = requests.get(url)
        if response.status_code == 200:
            return response.text.strip()
    except:
        return None
    return None

# Generate 3D molecule and export as .ply
def generate_ply(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    radii = {1: 0.3, 6: 0.5, 7: 0.45, 8: 0.4}
    colors = {1: "white", 6: "gray", 7: "blue", 8: "red"}

    plotter = pv.Plotter(off_screen=True)
    points = []
    labels = []
    mesh = None

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        x, y, z = pos.x, pos.y, pos.z
        rad = radii.get(atom.GetAtomicNum(), 0.4)
        color = colors.get(atom.GetAtomicNum(), "green")
        sphere = pv.Sphere(radius=rad, center=[x, y, z])
        plotter.add_mesh(sphere, color=color)
        plotter.add_point_labels([[x, y, z]], [atom.GetSymbol()], font_size=12, text_color="black")

    # Export the 3D mesh as .ply
    with tempfile.NamedTemporaryFile(delete=False, suffix=".ply") as tmp:
        plotter.screenshot()  # Forces mesh computation
        meshes = plotter.actors
        all_mesh = None
        for actor in meshes.values():
            if actor.mapper and hasattr(actor.mapper, 'dataset'):
                mesh_data = actor.mapper.dataset
                if all_mesh is None:
                    all_mesh = mesh_data
                else:
                    all_mesh += mesh_data
        if all_mesh:
            all_mesh.save(tmp.name)
            return tmp.name
    return None

# Input
user_input = st.text_input("🔹 Molecule name or SMILES:")

if user_input:
    if any(char in user_input for char in "=#[]()123456789"):
        smiles = user_input
    else:
        smiles = get_smiles_from_name(user_input)

    if smiles:
        st.success(f"Found SMILES: {smiles}")
        if st.button("🔄 Generate .ply File"):
            ply_path = generate_ply(smiles)
            if ply_path:
                with open(ply_path, "rb") as file:
                    st.download_button(
                        label="📥 Download .ply File",
                        data=file,
                        file_name="molecule.ply",
                        mime="application/octet-stream"
                    )
            else:
                st.error("⚠️ Error saving .ply file.")
    else:
        st.error("❌ SMILES not found or invalid input.")
