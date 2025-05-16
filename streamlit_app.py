import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import requests

# Function to get SMILES from molecule name using PubChem
def get_smiles(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        response = requests.get(url)
        smiles = response.json()['PropertyTable']['Properties'][0]['IsomericSMILES']
        return smiles
    except Exception:
        return None

# Visualize molecule in 3D with colors and labels
def visualize_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    AllChem.UFFOptimizeMolecule(mol)

    mol_block = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=500, height=400)
    view.addModel(mol_block, 'mol')
    # Show ball and stick style with atom colors
    view.setStyle({'stick': {}, 'sphere': {'radius': 0.3}})
    view.setBackgroundColor('0xeeeeee')
    view.zoomTo()

    # Add atom labels
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

        # Approximate lone pairs for N, O, S, Cl (just visual)
        if atom.GetAtomicNum() in [7, 8, 16, 17]:
            view.addLabel(
                "•",
                {'position': {'x': pos.x + 0.4, 'y': pos.y + 0.4, 'z': pos.z + 0.4},
                 'backgroundColor': 'lightblue',
                 'fontSize': 16,
                 'fontColor': 'blue'}
            )

    return view

st.title("Molecule Visualizer with Py3Dmol")

user_input = st.text_input("Enter molecule name or SMILES")

if user_input:
    # Decide if input is SMILES or name
    if all(c in 'CNOPSHFIBrcl1234567890-=#()@+[]' for c in user_input):  
        smiles = user_input
    else:
        smiles = get_smiles(user_input)

    if smiles:
        st.success(f"SMILES: {smiles}")
        view = visualize_molecule(smiles)
        st.components.v1.html(view._make_html(), height=450)
    else:
        st.error("Could not find SMILES for the input. Please try again.")
