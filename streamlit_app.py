import streamlit as st
import py3Dmol
import requests

st.set_page_config(page_title="3D Molecule Visualizer", layout="wide")

st.title("🔬 3D Molecule Visualizer (Name or SMILES)")
st.write("Enter a molecule name (e.g., water, glucose) or a SMILES string (e.g., CCO)")

user_input = st.text_input("Molecule Name or SMILES", "")

def fetch_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/TXT"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    return None

def show_molecule(smiles):
    mol = py3Dmol.view(width=500, height=500)
    mol.addModel(f'{smiles}', 'smi')
    mol.setStyle({'stick': {}})
    mol.zoomTo()
    mol.show()
    return mol

if user_input:
    if any(c in user_input for c in ['=', '#', '(', ')', '[', ']', '@', '+', '-', '/', '\\']) or len(user_input) > 20:
        # Looks like a SMILES
        smiles = user_input
    else:
        # Try converting name to SMILES
        smiles = fetch_smiles_from_name(user_input)
        if not smiles:
            st.error("Could not find a molecule for that name.")
        else:
            st.success(f"Found SMILES: {smiles}")

    if smiles:
        st.subheader("🧪 3D Structure")
        mol = show_molecule(smiles)
        mol.show()
        py3Dmol.render(mol)

# VR/AR section
       st.markdown("### 🌐 Want to see it in AR/VR?")
       ar_url = f"https://3dviewer.net/#modelurl=https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/SDF"
        st.markdown(f"[Click here to view in 3D AR/VR Viewer 🚀]({ar_url})")
