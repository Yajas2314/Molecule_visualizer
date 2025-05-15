import streamlit as st
import requests
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

st.set_page_config(layout="wide")
st.title("🔬 Enhanced 3D Molecule Visualizer")

# Input type
input_type = st.radio("Input Type", ("Molecule Name", "SMILES"))
user_input = st.text_input(f"Enter the {input_type}:")

# Get SMILES from molecule name
def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/TXT"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    return None

# Show 3D Molecule with full atoms and lone pairs
def show_3d_molecule_with_lone_pairs(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add explicit hydrogens
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    # Get XYZ coordinates
    conf = mol.GetConformer()
    xyz = ""
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"

        # Simulate lone pairs (O, N, F, Cl)
        if atom.GetSymbol() in ["O", "N", "F", "Cl"]:
            # Add 1–2 dots near these atoms
            lp_x = pos.x + 0.3
            lp_y = pos.y + 0.3
            lp_z = pos.z
            xyz += f"X {lp_x:.4f} {lp_y:.4f} {lp_z:.4f}\n"

    # Visualize in py3Dmol
    view = py3Dmol.view(width=800, height=500)
    view.addModel(xyz, "xyz")
    view.setStyle({
        "stick": {},
        "sphere": {"scale": 0.3},
        "nonbonded": {"color": "gray"}
    })
    view.setBackgroundColor("white")
    view.zoomTo()
    return view

# MAIN APP
if user_input:
    if input_type == "Molecule Name":
        smiles = get_smiles_from_name(user_input)
        if not smiles:
            st.error("Could not find SMILES for this molecule.")
        else:
            st.success(f"SMILES: {smiles}")
    else:
        smiles = user_input

    if smiles:
        st.subheader("🧬 3D Structure with Lone Pairs")
        viewer = show_3d_molecule_with_lone_pairs(smiles)
        html = viewer._make_html()
        st.components.v1.html(html, height=500, width=800)
        st.info("Atoms, Hydrogens, and Lone Pair dots are visualized.")
