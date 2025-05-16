import streamlit as st
import py3Dmol
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

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

def draw_3d_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES input. Please enter a correct SMILES string or molecule name.")
    
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # Get 3D coordinates
    mol_block = Chem.MolToMolBlock(mol)
    
    view = py3Dmol.view(width=700, height=500)
    view.addModel(mol_block, "mol")
    
    # Stick & Ball
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    
    # Element labels
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        element = atom.GetSymbol()
        view.addLabel(
            element,
            {'position': {'x': pos.x, 'y': pos.y, 'z': pos.z},
             'backgroundColor': 'white', 'fontColor': 'black', 'fontSize': 12}
        )

        # Lone pair visualization (approximate)
        if atom.GetAtomicNum() in [7, 8, 9, 16]:  # N, O, F, S
            lp_offset = 0.5
            view.addSphere({
                'center': {'x': pos.x + lp_offset, 'y': pos.y, 'z': pos.z},
                'radius': 0.1,
                'color': 'gray',
                'opacity': 0.6
            })

    view.zoomTo()
    return view

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
        viewer = draw_3d_molecule(smiles)
        st.components.v1.html(viewer._make_html(), height=500)
