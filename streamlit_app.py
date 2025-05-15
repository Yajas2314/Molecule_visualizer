import streamlit as st
import py3Dmol
import requests

# --- Page settings ---
st.set_page_config(page_title="3D Molecule Viewer", layout="centered")
st.title("🧪 3D Molecule Visualizer")

# --- Function to get SMILES from molecule name ---
def get_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        response = requests.get(url)
        data = response.json()
        return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except Exception:
        return None

# --- Function to generate 3D molecule viewer HTML ---
def generate_3d_view(smiles):
    view = py3Dmol.view(width=500, height=500)
    view.addModel(smiles, 'smi')
    view.setStyle({'stick': {}})
    view.zoomTo()
    return view._make_html()

# --- Sidebar input ---
st.sidebar.header("Input Molecule")
molecule_name = st.sidebar.text_input("Enter Molecule Name (e.g. glucose):")
smiles_input = st.sidebar.text_input("Or Enter SMILES Code (e.g. C(CO)O):")

# --- Choose input ---
if smiles_input:
    smiles = smiles_input.strip()
elif molecule_name:
    smiles = get_smiles_from_name(molecule_name)
else:
    smiles = None

# --- Display ---
if smiles:
    st.subheader("🔬 3D Structure Viewer")
    html = generate_3d_view(smiles)
    with st.container():
        st.components.v1.html(html, height=520)
    st.code(smiles, language='text')
else:
    st.info("Please enter a molecule name or SMILES code in the sidebar.")
