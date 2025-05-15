import streamlit as st
import py3Dmol
import requests

# --- Streamlit page setup ---
st.set_page_config(page_title="3D Molecule Viewer", layout="wide")
st.title("🧪 3D Molecule Visualizer")
st.subheader("🧬 3D Structure Viewer")

# --- SMILES fetcher from name ---
def get_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        res = requests.get(url)
        res.raise_for_status()
        return res.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except Exception:
        return None

# --- 3D Viewer from SMILES ---
def generate_3d_view(smiles):
    view = py3Dmol.view(width=500, height=500)
    view.addModel(smiles, 'smi')
    view.setStyle({'stick': {}})
    view.zoomTo()
    return view._make_html()

# --- Sidebar inputs ---
st.sidebar.header("Input Molecule")
mol_name = st.sidebar.text_input("Enter Molecule Name (e.g. glucose):")
smiles_input = st.sidebar.text_input("Or Enter SMILES Code (e.g. C(CO)O):")

# --- Get SMILES code ---
if smiles_input:
    smiles = smiles_input.strip()
elif mol_name:
    smiles = get_smiles_from_name(mol_name.strip())
else:
    smiles = None

# --- Display area directly below heading ---
if smiles:
    col1, col2, col3 = st.columns([1, 2, 1])  # Centered layout
    with col2:
        html = generate_3d_view(smiles)
        st.components.v1.html(html, height=520)
    st.success(f"SMILES Code: `{smiles}`")
else:
    st.info("Please enter a molecule name or SMILES code in the sidebar.")
