import streamlit as st
import py3Dmol
from urllib.parse import quote
import requests

# Get SMILES code from a molecule name using PubChem API
def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(name)}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        try:
            data = response.json()
            smiles = data['PropertyTable']['Properties'][0]['IsomericSMILES']
            return smiles
        except Exception:
            return None
    return None

# Show molecule with enhanced styles (ball-and-stick, colors)
def show_molecule(smiles):
    view = py3Dmol.view(width=500, height=500)
    view.addModel(smiles, "smi")
    # Set enhanced style: ball-and-stick, colored by atom
    style = {
        "stick": {"radius": 0.2},
        "sphere": {"scale": 0.3},
        "cartoon": {},
        "atom": {"colorscheme": "Jmol"}
    }
    view.setStyle({}, style)
    view.zoomTo()
    return view

# Streamlit App Interface
st.title("🧪 Advanced 3D Molecule Visualizer")
st.markdown("""
This app lets you visualize any molecule in **3D** with **enhanced styles**.
- Input a molecule name or a SMILES code.
- The structure is rendered in ball-and-stick format.
- Colors distinguish different atoms.
- Zoom, rotate, and inspect in 3D.
""")

# Input method selection
input_type = st.radio("Choose input method:", ("Molecule Name", "SMILES Code"))

smiles = ""

if input_type == "Molecule Name":
    name = st.text_input("Enter molecule name:")
    if name:
        smiles = get_smiles_from_name(name)
        if smiles:
            st.success(f"Found SMILES: `{smiles}`")
        else:
            st.error("Could not fetch SMILES for the given molecule name.")
else:
    smiles = st.text_input("Enter SMILES code:")

# Render molecule if SMILES is available
if smiles:
    st.subheader("🧬 3D Structure")
    mol = show_molecule(smiles)
    mol_html = mol._make_html()
    st.components.v1.html(mol_html, height=500)

    st.info("Use your mouse to zoom, rotate, and explore the structure.")

# Note about lone pairs
st.markdown("""
---
**Note**: Visualization of *lone pair electrons* is not directly supported by 3Dmol.js.
However, you can infer electron domains from the structure geometry.
""")
