import streamlit as st
import py3Dmol
import requests

# --- Page settings ---
st.set_page_config(page_title="3D Molecule Viewer", layout="wide")
st.title("🧪 3D Molecule Visualizer")

# --- Function to get SMILES from molecule name using PubChem ---
def get_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        response = requests.get(url)
        data = response.json()
        return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except Exception:
        return None

# --- Function to generate 3D HTML view from SMILES ---
def generate_3d_view(smiles):
    view = py3Dmol.view(width=500, height=500)
    view.addModel(smiles, 'smi')
    view.setStyle({'stick': {}})
    view.zoomTo()
    return view._make_html()

# --- Sidebar inputs ---
st.sidebar.header("Input")
molecule_name = st.sidebar.text_input("Enter Molecule Name (e.g. glucose):")
smiles_input = st.sidebar.text_input("Or Enter SMILES Code (e.g. C(CO)O):")

# --- Determine which input to use ---
if smiles_input:
    smiles = smiles_input.strip()
elif molecule_name:
    smiles = get_smiles_from_name(molecule_name)
else:
    smiles = None

# --- Output section ---
if smiles:
    st.subheader("🔬 3D Structure")
    try:
        mol_html = generate_3d_view(smiles)
        st.components.v1.html(mol_html, height=520)
        st.success(f"SMILES: {smiles}")
    except Exception as e:
        st.error(f"Could not render molecule. Error: {e}")
else:
    st.info("Please enter a molecule name or SMILES code in the sidebar.")
