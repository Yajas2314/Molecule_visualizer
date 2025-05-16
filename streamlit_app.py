import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import pyvista as pv
import tempfile

def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/TXT"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    return None

def generate_ply(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    radii = {
        1: 0.25, 6: 0.70, 7: 0.65, 8: 0.60,
        9: 0.50, 15: 1.00, 16: 1.00, 17: 0.85,
    }

    spheres = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        radius = radii.get(atom.GetAtomicNum(), 0.5)
        sphere = pv.Sphere(radius=radius, center=[pos.x, pos.y, pos.z])
        spheres.append(sphere)

    combined = spheres[0]
    for sphere in spheres[1:]:
        combined += sphere

    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".ply")
    combined.save(temp_file.name)
    return temp_file.name

st.title("3D Molecule Visualizer with AR (.glb)")
user_input = st.text_input("Enter SMILES or Molecule Name")

if user_input:
    if any(char.isalpha() for char in user_input) and not any(char in user_input for char in "#=()[]"):
        smiles = get_smiles_from_name(user_input)
    else:
        smiles = user_input

    if smiles:
        st.success(f"Found SMILES: {smiles}")
        if st.button("Generate .ply File"):
            ply_path = generate_ply(smiles)
            with open(ply_path, "rb") as f:
                st.download_button("Download .ply", f, file_name="molecule.ply")
            
    else:
        st.error("Invalid molecule name or SMILES")
