import streamlit as st
import py3Dmol
import requests

st.set_page_config(page_title="3D Molecule Visualizer", layout="wide")
st.title("🧪 3D Molecule Visualizer")

# Sidebar inputs
st.sidebar.header("Input Molecule")
mol_name = st.sidebar.text_input("Enter Molecule Name (e.g. glucose):")
smiles_input = st.sidebar.text_input("Or Enter SMILES Code (e.g. C(CO)O):")
style = st.sidebar.selectbox("Choose Visualization Style:", ["Ball and Stick", "Stick", "Sphere", "Line"])

def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        return data['PropertyTable']['Properties'][0]['CanonicalSMILES']
    except Exception as e:
        st.error(f"Could not fetch SMILES for molecule name '{name}': {e}")
        return None

def show_molecule(smiles, style):
    style_dict = {
        "Ball and Stick": {'stick': {}, 'sphere': {'scale': 0.3}},
        "Stick": {'stick': {}},
        "Sphere": {'sphere': {}},
        "Line": {'line': {}}
    }
    view = py3Dmol.view(width=700, height=500)
    view.addModel(smiles, 'smi')
    view.addHydrogens()
    view.setStyle(style_dict[style])
    view.setBackgroundColor('0xeeeeee')
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=520)

# Determine SMILES code from input
smiles = None
if smiles_input.strip():
    smiles = smiles_input.strip()
elif mol_name.strip():
    smiles = get_smiles_from_name(mol_name.strip())

if smiles:
    st.markdown("### 3D Structure Viewer:")
    show_molecule(smiles, style)
    st.success(f"Showing structure for: `{smiles}`")
else:
    st.info("Enter a molecule name or SMILES code to display its 3D structure.")
