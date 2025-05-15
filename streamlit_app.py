import streamlit as st
import py3Dmol
import requests

# Page configuration
st.set_page_config(page_title="3D Molecule Visualizer", layout="wide")

# Title and subtitle
st.title("🧪 3D Molecule Visualizer")
st.markdown("### 🔬 3D Structure Viewer")

# Sidebar input
st.sidebar.header("Input Molecule")
mol_name = st.sidebar.text_input("Enter Molecule Name (e.g. glucose)")
smiles_input = st.sidebar.text_input("Or Enter SMILES Code (e.g. C(CO)O)")

# Style selector
style = st.sidebar.selectbox("Choose Visualization Style:", ["Ball and Stick", "Stick", "Sphere", "Line"])

# Function to fetch SMILES from name
def fetch_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        response = requests.get(url)
        data = response.json()
        return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except:
        return None

# Function to generate 3D HTML from SMILES
def generate_molecule_html(smiles_code, style_type):
    view = py3Dmol.view(width=700, height=500)
    view.addModel(smiles_code, "smi")
    view.addHydrogens()
    view.setBackgroundColor("white")

    if style_type == "Ball and Stick":
        view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    elif style_type == "Stick":
        view.setStyle({'stick': {}})
    elif style_type == "Sphere":
        view.setStyle({'sphere': {}})
    elif style_type == "Line":
        view.setStyle({'line': {}})

    view.zoomTo()
    return view._make_html()

# Determine which SMILES to use
if smiles_input:
    smiles = smiles_input.strip()
elif mol_name:
    smiles = fetch_smiles_from_name(mol_name.strip())
else:
    smiles = None

# MAIN AREA: show molecule if available
if smiles:
    # Show the 3D viewer directly under heading
    st.markdown("#### 🧬 3D Structure Below:")
    html = generate_molecule_html(smiles, style)
    st.components.v1.html(html, height=520)
    st.success(f"**SMILES used:** `{smiles}`")
else:
    st.info("🔍 Please enter a molecule name or SMILES code in the sidebar to visualize.")
