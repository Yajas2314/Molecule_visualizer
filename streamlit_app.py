import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import requests
import tempfile
import pyvista as pv

# Get SMILES from compound name via PubChem
def get_smiles(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        response = requests.get(url)
        smiles = response.json()['PropertyTable']['Properties'][0]['IsomericSMILES']
        return smiles
    except Exception:
        return None

# Visualize molecule with py3Dmol
def visualize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    AllChem.UFFOptimizeMolecule(mol)

    mol_block = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=500, height=400)
    view.addModel(mol_block, 'mol')
    view.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
    view.setBackgroundColor('0xeeeeee')
    view.zoomTo()

    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        view.addLabel(
            atom.GetSymbol(),
            {'position': {'x': pos.x, 'y': pos.y, 'z': pos.z},
             'backgroundColor': 'white',
             'fontSize': 14,
             'fontColor': 'black',
             'borderThickness': 1,
             'borderColor': 'gray'}
        )

        if atom.GetAtomicNum() in [7, 8, 16, 17]:
            view.addLabel(
                "•",
                {'position': {'x': pos.x + 0.4, 'y': pos.y + 0.4, 'z': pos.z + 0.4},
                 'backgroundColor': 'lightblue',
                 'fontSize': 16,
                 'fontColor': 'blue'}
            )
    return view

import numpy as np

def generate_glb(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    # Atomic radii (approximate, in angstroms)
    atomic_radii = {
        1: 0.25,  # H
        6: 0.70,  # C
        7: 0.65,  # N
        8: 0.60,  # O
        9: 0.50,  # F
        15: 1.00, # P
        16: 1.00, # S
        17: 0.85, # Cl
        # Add more if needed
    }

    # Create pyvista sphere meshes for each atom
    spheres = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        radius = atomic_radii.get(atom.GetAtomicNum(), 0.5)
        sphere = pv.Sphere(radius=radius, center=[pos.x, pos.y, pos.z], theta_resolution=16, phi_resolution=16)
        spheres.append(sphere)

    # Combine all spheres into a single mesh
    combined = spheres[0]
    for s in spheres[1:]:
        combined = combined + s

    # Save to temporary .glb file
    tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".glb")
    try:
        combined.save(tmp_file.name)
        return tmp_file.name
    except Exception as e:
        st.error(f"Error saving .glb file: {e}")
        return None
    

# Streamlit app starts here
st.title("Molecule Visualizer with 3Dmol.js and .glb Export")

user_input = st.text_input("Enter molecule name or SMILES")

if user_input:
    # Check if input looks like SMILES or name
    if all(c in 'CNOPSHFIBrcl1234567890-=#()@+[]' for c in user_input):
        smiles = user_input
    else:
        smiles = get_smiles(user_input)

    if smiles:
        st.success(f"Using SMILES: {smiles}")

        # Show molecule 3D viewer
        view = visualize_molecule(smiles)
        st.components.v1.html(view._make_html(), height=450)

        # Button to generate .glb file
        if st.button("Generate .glb file"):
            glb_path = generate_glb(smiles)
            if glb_path:
                with open(glb_path, "rb") as f:
                    glb_bytes = f.read()
                st.download_button(
                    label="Download .glb 3D model",
                    data=glb_bytes,
                    file_name="molecule_model.glb",
                    mime="model/gltf-binary"
                )
    else:
        st.error("Could not find SMILES for the input. Please try again.")
