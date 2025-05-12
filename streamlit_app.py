import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import requests

st.set_page_config(page_title="3D Molecule Visualizer with VR", layout="centered")
st.title("🔬 3D Molecule Visualizer with VR Support")

# Input section
option = st.radio("Choose input type:", ("Molecule Name", "SMILES Code"))

if option == "Molecule Name":
    molecule_name = st.text_input("Enter Molecule Name (e.g. ethanol):")
    if molecule_name:
        # Fetch SMILES using PubChem API
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{molecule_name}/property/IsomericSMILES/JSON"
            response = requests.get(url)
            data = response.json()
            smiles = data['PropertyTable']['Properties'][0]['IsomericSMILES']
            st.success(f"Found SMILES: {smiles}")
        except Exception as e:
            st.error("Could not find SMILES for the given name. Please check spelling.")
            st.stop()
else:
    smiles = st.text_input("Enter SMILES Code (e.g. CCO):")

if smiles:
    st.subheader("🧪 3D Structure")

    def show_molecule(smiles_code):
        mol = py3Dmol.view(width=500, height=400)
        mol.addModel(smiles_code, 'smi')
        mol.setStyle({"stick": {}})
        mol.zoomTo()
        return mol

    mol = show_molecule(smiles)
    components.html(mol._make_html(), height=400)

    st.markdown("---")
    st.subheader("🌐 VR/AR View (Experimental)")
    st.markdown("You can also explore the 3D structure in *augmented/virtual reality* via the site below.")
    ar_vr_url = f"https://3dviewer.net/?load=https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/record/SDF/?record_type=3d"
    st.markdown(f"[Open in AR/VR Viewer 🔗]({ar_vr_url})")
