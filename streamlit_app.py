import streamlit as st
import py3Dmol
import requests

st.set_page_config(page_title="3D Molecule Viewer", layout="wide")

# Header
st.title("🧪 3D Molecule Visualizer")
st.markdown("### 🔬 3D Structure Viewer")

# Sidebar input
st.sidebar.header("Input Molecule")
mol_name = st.sidebar.text_input("Enter Molecule Name (e.g. glucose)")
smiles_input = st.sidebar.text_input("Or Enter SMILES Code (e.g. C(CO)O)")

# Style selector
style = st.sidebar.selectbox("Choose Visualization Style:", ["Ball and Stick", "Stick", "Sphere", "Line"])

# Function: Fetch SMILES from name
def fetch_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        response = requests.get(url)
        data = response.json()
        return data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
    except:
        return None

# Function: Create 3D HTML
def generate_molecule_html(smiles_code, style_type):
    view = py3Dmol.view(width=600, height=500)
    view.addModel(smiles_code, "smi")
    view.addHydrogens()  # Show explicit H
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

# SMILES decision
if smiles_input:
    smiles = smiles_input.strip()
elif mol_name:
    smiles = fetch_smiles_from_name(mol_name.strip())
else:
    smiles = None

# Show 3D viewer
if smiles:
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.markdown("#### 🧪 Molecule Structure:")
        html = generate_molecule_html(smiles, style)
        st.components.v1.html(html, height=520)
        st.success(f"**SMILES Used:** `{smiles}`")
else:
    st.info("Please enter a molecule name or SMILES code in the sidebar.")
