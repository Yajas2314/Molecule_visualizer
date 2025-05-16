import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import requests
import trimesh

# Get SMILES from compound name
def get_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        response = requests.get(url)
        smiles = response.json()['PropertyTable']['Properties'][0]['IsomericSMILES']
        return smiles
    except:
        return None

# Show molecule in 3D
def show_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    view = py3Dmol.view(width=500, height=400)
    view.addModel(Chem.MolToMolBlock(mol), 'mol')
    view.setStyle({'stick': {}})
    view.setBackgroundColor('white')
    view.zoomTo()

    # Atom labels and approximate lone pair marker
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        view.addLabel(atom.GetSymbol(), {
            'position': {'x': pos.x, 'y': pos.y, 'z': pos.z},
            'backgroundColor': 'white',
            'fontSize': 14
        })
        # Add "•" as lone pair approximation (visual only, not accurate)
        if atom.GetAtomicNum() in [7, 8, 16, 17]:  # N, O, S, Cl
            view.addLabel("•", {
                'position': {'x': pos.x + 0.5, 'y': pos.y + 0.5, 'z': pos.z + 0.5},
                'backgroundColor': 'lightblue',
                'fontSize': 12
            })

    return view

# Generate .glb file
def generate_glb_from_smiles(smiles, filename):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    vertices = []
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        vertices.append([pos.x, pos.y, pos.z])

    mesh = trimesh.points.PointCloud(vertices)
    mesh.export(filename)
    return True

# Streamlit UI
st.title("🧪 Molecule Visualizer with AR (.glb)")

user_input = st.text_input("Enter Molecule Name or SMILES Code")

if user_input:
    if any(c in user_input for c in "=#()"):  # likely SMILES
        smiles = user_input
    else:
        smiles = get_smiles_from_name(user_input)

    if smiles:
        st.success(f"Found SMILES: {smiles}")

        st.subheader("🔬 3D Structure Viewer")
        view = show_molecule(smiles)
        st.components.v1.html(view._make_html(), height=450)

        st.subheader("📦 Export .glb for AR")
        if st.button("Generate .glb File"):
            path = "molecule_model.glb"
            if generate_glb_from_smiles(smiles, path):
                with open(path, "rb") as f:
                    st.download_button("Download .glb File", f, "molecule_model.glb")
            else:
                st.error("Could not generate the .glb file.")
    else:
        st.error("SMILES not found. Try a different name or valid SMILES string.")
