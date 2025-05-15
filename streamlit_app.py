import streamlit as st
import py3Dmol
import requests

st.set_page_config(page_title="3D Molecule Visualizer", layout="wide")

st.title("🧪 3D Molecule Visualizer")
st.markdown("### 🔬 3D Structure Viewer")

# Sidebar Inputs
st.sidebar.header("Input Molecule")
mol_name = st.sidebar.text_input("Enter Molecule Name (e.g. glucose)")
smiles_input = st.sidebar.text_input("Or Enter SMILES Code (e.g. C(CO)O)")
style = st.sidebar.selectbox("Choose Visualization Style:", ["Ball and Stick", "Stick", "Sphere", "Line"])

# Function to fetch SMILES from name
def get_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
        res = requests.get(url)
        res.raise_for_status()
        data = res.json()
        return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except:
        return None

# Function to generate 3Dmol HTML
def show_molecule(smiles, style):
    view = py3Dmol.view(width=700, height=500)
    view.addModel(smiles, 'smi')
    view.addHydrogens()
    view.setBackgroundColor("white")

    # Apply style
    if style == "Ball and Stick":
        view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    elif style == "Stick":
        view.setStyle({'stick': {}})
    elif style == "Sphere":
        view.setStyle({'sphere': {}})
    elif style == "Line":
        view.setStyle({'line': {}})
    
    view.zoomTo()
    return view._make_html()

# Determine SMILES
smiles = ""
if smiles_input:
    smiles = smiles_input.strip()
elif mol_name:
    smiles = get_smiles_from_name(mol_name.strip())

# Render 3D Structure
if smiles:
    st.markdown("### 🧬 3D Structure Below:")
    html = show_molecule(smiles, style)
    st.components.v1.html(html, height=520)
    st.success(f"✅ Molecule Loaded: `{smiles}`")
else:
    st.info("ℹ️ Please enter a molecule name or SMILES code to display its 3D structure.")
