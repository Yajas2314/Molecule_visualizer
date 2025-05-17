import streamlit as st
import py3Dmol
import requests


st.set_page_config(page_title="Molecule Viewer", layout="wide")
st.title("🔬 Molecule Visualizer")

# Get user input
user_input = st.text_input("Enter SMILES or Molecule Name", "")

def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    try:
        response = requests.get(url)
        smiles = response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
        return smiles
    except:
        return None
from rdkit import Chem
smiles = "CCO"
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

from rdkit.Chem import AllChem
AllChem.EmbedMolecule(mol)

mol_block = Chem.MolToMolBlock(mol)

def draw_3d_molecule_with_lone_pairs(mol):
    mb = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(mb, 'mol')
    viewer.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    viewer.zoomTo()

    # Add approximate lone pair spheres near O or N atoms
    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['O', 'N']:
            pos = conf.GetAtomPosition(atom.GetIdx())
            viewer.addSphere({
                'center': {'x': pos.x + 0.3, 'y': pos.y + 0.3, 'z': pos.z},
                'radius': 0.15,
                'color': 'gray',
                'opacity': 0.6
            })
            viewer.addSphere({
                'center': {'x': pos.x - 0.3, 'y': pos.y - 0.3, 'z': pos.z},
                'radius': 0.15,
                'color': 'gray',
                'opacity': 0.6
            })

    return viewer


# Logic to convert molecule name to SMILES
if user_input:
    if all(c.isalpha() or c.isdigit() for c in user_input):
        smiles = get_smiles_from_name(user_input)
        if not smiles:
            st.error("Could not fetch SMILES. Please check the molecule name.")
        else:
            st.success(f"Found SMILES: {smiles}")
    else:
        smiles = user_input

    if smiles:
        viewer = draw_3d_molecule_with_lone_pairs(smiles)
        st.components.v1.html(viewer._make_html(), height=500)

ply_path = generate_colored_labeled_ply("CCO")  # Ethanol
print("Saved to:", ply_path)
