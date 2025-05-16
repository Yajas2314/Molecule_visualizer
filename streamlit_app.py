import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import requests
import pyvista as pv
import os

st.set_page_config(layout="wide")
st.title("🔬 3D Molecule Viewer with Atom Colors, Labels & Lone Pairs")

input_type = st.radio("Input Type", ("Molecule Name", "SMILES"))
user_input = st.text_input(f"Enter the {input_type}:")

# Atom-specific colors
atom_colors = {
    "H": "white",
    "C": "gray",
    "N": "blue",
    "O": "red",
    "F": "green",
    "Cl": "green",
    "Br": "brown",
    "I": "purple",
    "S": "yellow",
    "P": "orange"
}

def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/TXT"
    res = requests.get(url)
    if res.status_code == 200:
        return res.text.strip()
    return None

def generate_3d_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    view = py3Dmol.view(width=800, height=500)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        sym = atom.GetSymbol()
        color = atom_colors.get(sym, "lightgray")

        # Atom as sphere
        view.addSphere({
            'center': {'x': pos.x, 'y': pos.y, 'z': pos.z},
            'radius': 0.3,
            'color': color,
            'opacity': 1.0
        })

        # Atom label
        view.addLabel(sym, {
            'position': {'x': pos.x, 'y': pos.y, 'z': pos.z},
            'backgroundColor': 'black',
            'fontColor': 'white',
            'fontSize': 10
        })

        # Lone pairs for O, N, F, Cl
        if sym in ["O", "N", "F", "Cl"]:
            lp1 = (pos.x + 0.4, pos.y + 0.3, pos.z)
            lp2 = (pos.x - 0.3, pos.y - 0.3, pos.z)
            view.addSphere({
                'center': {'x': lp1[0], 'y': lp1[1], 'z': lp1[2]},
                'radius': 0.1,
                'color': 'lightgray',
                'opacity': 1.0
            })
            view.addSphere({
                'center': {'x': lp2[0], 'y': lp2[1], 'z': lp2[2]},
                'radius': 0.1,
                'color': 'lightgray',
                'opacity': 1.0
            })

    # Bonds as sticks
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        begin_pos = conf.GetAtomPosition(begin_idx)
        end_pos = conf.GetAtomPosition(end_idx)
        view.addCylinder({
            'start': {'x': begin_pos.x, 'y': begin_pos.y, 'z': begin_pos.z},
            'end': {'x': end_pos.x, 'y': end_pos.y, 'z': end_pos.z},
            'radius': 0.1,
            'color': 'gray',
            'opacity': 1.0
        })

    view.setBackgroundColor('white')
    view.zoomTo()
    return view

if user_input:
    if input_type == "Molecule Name":
        smiles = get_smiles_from_name(user_input)
        if not smiles:
            st.error("❌ Could not find SMILES for this molecule.")
        else:
            st.success(f"✅ SMILES: {smiles}")
    else:
        smiles = user_input.strip()

    if smiles:
        st.subheader("🧬 3D Structure")
        mol_view = generate_3d_structure(smiles)
        mol_html = mol_view._make_html()
        st.components.v1.html(mol_html, height=500, width=800)

def generate_glb_from_smiles(smiles: str, output_file: str = "molecule.glb"):
    # Generate RDKit molecule with 3D coordinates
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    plotter = pv.Plotter(off_screen=True)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        sphere = pv.Sphere(radius=0.3, center=(pos.x, pos.y, pos.z))
        plotter.add_mesh(sphere, color='white')

    mesh = plotter.mesh
    if mesh is None:
        print("Error: no mesh generated")
        return

if st.button("Generate AR Model"):
    try:
        generate_glb_from_smiles(smiles)
        with open("molecule.glb", "rb") as f :
            st.download_button("Download .glb file", f, "molecule.glb")
        st.markdown("[Click here to view in AR](https://your-webar-viewer.com/molecule.glb)")
    except Exception as e:
        st.error(f"Error generating model:{e}")

    # Save to .glb file
    mesh.save(output_file)
    print(f"Saved molecule model as {output_file}")
