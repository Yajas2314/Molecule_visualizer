import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

# Function to fetch SMILES from molecule name using PubChem
def fetch_smiles(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/TXT"
        response = requests.get(url)
        if response.status_code == 200:
            return response.text.strip()
    except:
        return None

# Function to create 3D molecule with labels and lone pair dots
def draw_molecule_with_labels_and_lone_pairs(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    mol_block = Chem.MolToMolBlock(mol)

    viewer = py3Dmol.view(width=600, height=500)
    viewer.addModel(mol_block, 'mol')

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = mol.GetConformer().GetAtomPosition(idx)
        label = atom.GetSymbol()
        color = {
            'H': 'white', 'C': 'gray', 'O': 'red', 'N': 'blue',
            'S': 'yellow', 'P': 'orange', 'F': 'green', 'Cl': 'green',
            'Br': 'brown', 'I': 'purple'
        }.get(label, 'black')

        viewer.addSphere({'center': {'x': pos.x, 'y': pos.y, 'z': pos.z},
                          'radius': 0.3, 'color': color})
        viewer.addLabel(label, {
            'position': {'x': pos.x, 'y': pos.y, 'z': pos.z},
            'backgroundColor': 'black', 'fontColor': 'white', 'fontSize': 12
        })

        # Add approximate lone pairs
        lone_pairs = {
            'O': 2, 'N': 1, 'S': 2, 'F': 3, 'Cl': 3, 'Br': 3, 'I': 3
        }.get(label, 0)
        for i in range(lone_pairs):
            offset = (i + 1) * 0.25
            viewer.addSphere({
                'center': {'x': pos.x + offset, 'y': pos.y + offset, 'z': pos.z + offset},
                'radius': 0.1, 'color': 'black'
            })

    viewer.setStyle({'stick': {}})
    viewer.setBackgroundColor('white')
    viewer.zoomTo()
    return viewer

# Streamlit App
st.title("Universal Molecule Visualizer")
user_input = st.text_input("Enter molecule name or SMILES:", "ethanol")

if user_input:
    if all(char.isalpha() or char in '=#()[]@+-' for char in user_input):
        smiles = fetch_smiles(user_input)
    else:
        smiles = user_input
    
    if smiles:
        st.success(f"SMILES: {smiles}")
        viewer = draw_molecule_with_labels_and_lone_pairs(smiles)
        html = viewer._make_html()
        st.components.v1.html(html, height=600, width=700)
    else:
        st.error("Could not fetch SMILES. Please try a valid molecule.")
