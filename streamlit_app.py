import streamlit as st
import py3Dmol
import requests

st.set_page_config(page_title="3D Molecule Visualizer", layout="wide")

st.title("🧪 3D Molecule Visualizer")
st.markdown("### 🔬 3D Structure Viewer")

# Sidebar inputs
st.sidebar.header("Input Molecule")
mol_name = st.sidebar.text_input("Enter Molecule Name (e.g. glucose)")
smiles_input = st.sidebar.text_input("Or Enter SMILES Code (e.g. C(CO)O)")
style = st.sidebar.selectbox("Choose Visualization Style:", ["Ball and Stick", "Stick", "Sphere", "Line"])

# Function to fetch SMILES from molecule name using PubChem API
def get_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        res = requests.get(url)
        res.raise_for_status()
        data = res.json()
        return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except Exception as e:
        st.error(f"Error fetching SMILES for '{name}': {e}")
        return None

# Function to display molecule using py3Dmol
def display_molecule(smiles, style):
    style_dict = {
        "Ball and Stick": {'stick': {}, 'sphere': {'scale': 0.3}},
        "Stick": {'stick': {}},
        "Sphere": {'sphere': {}},
        "Line": {'line': {}}
    }
    
    view = py3Dmol.view(width=800, height=500)
    view.addModel(smiles, 'smi')
    view.addHydrogens()
    view.setStyle(style_dict[style])
    view.setBackgroundColor("white")
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=520, scrolling=False)

# Get the SMILES string
smiles = None
if smiles_input.strip():
    smiles = smiles_input.strip()
elif mol_name.strip():
    smiles = get_smiles_from_name(mol_name.strip())

# Display the molecule or show info message
if smiles:
    st.markdown("### 🧬 3D Structure Below:")
    display_molecule(smiles, style)
    st.success(f"✅ Molecule Loaded: `{smiles}`")
else:
    st.info("ℹ️ Enter a molecule name or SMILES code to display its 3D structure.")
