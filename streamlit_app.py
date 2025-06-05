import streamlit as st
import requests
import tempfile
import base64
import os

# ---------- Static atomic data for elements 1-30 ----------
atomic_data = {
    'H':  {'name': 'Hydrogen',  'number': 1,  'mass': 1.008,  'config': '1s1'},
    'He': {'name': 'Helium',    'number': 2,  'mass': 4.0026, 'config': '1s2'},
    'Li': {'name': 'Lithium',   'number': 3,  'mass': 6.94,   'config': '1s2 2s1'},
    'Be': {'name': 'Beryllium', 'number': 4,  'mass': 9.0122, 'config': '1s2 2s2'},
    'B':  {'name': 'Boron',     'number': 5,  'mass': 10.81,  'config': '1s2 2s2 2p1'},
    'C':  {'name': 'Carbon',    'number': 6,  'mass': 12.011, 'config': '1s2 2s2 2p2'},
    'N':  {'name': 'Nitrogen',  'number': 7,  'mass': 14.007, 'config': '1s2 2s2 2p3'},
    'O':  {'name': 'Oxygen',    'number': 8,  'mass': 15.999, 'config': '1s2 2s2 2p4'},
    'F':  {'name': 'Fluorine',  'number': 9,  'mass': 18.998, 'config': '1s2 2s2 2p5'},
    'Ne': {'name': 'Neon',      'number': 10, 'mass': 20.180, 'config': '1s2 2s2 2p6'},
    'Na': {'name': 'Sodium',    'number': 11, 'mass': 22.990, 'config': '1s2 2s2 2p6 3s1'},
    'Mg': {'name': 'Magnesium', 'number': 12, 'mass': 24.305, 'config': '1s2 2s2 2p6 3s2'},
    'Al': {'name': 'Aluminium', 'number': 13, 'mass': 26.982, 'config': '1s2 2s2 2p6 3s2 3p1'},
    'Si': {'name': 'Silicon',   'number': 14, 'mass': 28.085, 'config': '1s2 2s2 2p6 3s2 3p2'},
    'P':  {'name': 'Phosphorus','number': 15, 'mass': 30.974, 'config': '1s2 2s2 2p6 3s2 3p3'},
    'S':  {'name': 'Sulfur',    'number': 16, 'mass': 32.06,  'config': '1s2 2s2 2p6 3s2 3p4'},
    'Cl': {'name': 'Chlorine',  'number': 17, 'mass': 35.45,  'config': '1s2 2s2 2p6 3s2 3p5'},
    'Ar': {'name': 'Argon',     'number': 18, 'mass': 39.948, 'config': '1s2 2s2 2p6 3s2 3p6'},
    'K':  {'name': 'Potassium', 'number': 19, 'mass': 39.098, 'config': '1s2 2s2 2p6 3s2 3p6 4s1'},
    'Ca': {'name': 'Calcium',   'number': 20, 'mass': 40.078, 'config': '1s2 2s2 2p6 3s2 3p6 4s2'},
    # Add more as needed up to 30
}

# ---------- Utility Functions ----------
def fetch_smiles(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/TXT"
    r = requests.get(url)
    return r.text.strip() if r.status_code == 200 else None

def generate_fake_glb(smiles):
    # Placeholder for actual glb generation logic
    # For now just encode the SMILES string as a dummy glb
    glb_filename = f"molecule_{smiles}.glb"
    with open(glb_filename, 'wb') as f:
        f.write(smiles.encode('utf-8'))
    return glb_filename

def get_element_info(smiles):
    atoms = {}
    for symbol in atomic_data:
        if symbol in smiles:
            atoms[symbol] = atomic_data[symbol]
    return atoms

def render_model_viewer(glb_url):
    viewer = f'''
    <model-viewer src="{glb_url}" ar ar-modes="webxr scene-viewer quick-look" auto-rotate camera-controls style="width: 100%; height: 500px;"></model-viewer>
    '''
    return viewer

# ---------- Streamlit App ----------
st.set_page_config(page_title="AR Molecule Visualizer", layout="wide")
st.title("🔬 AR Molecule Visualizer")

input_text = st.text_input("Enter molecule name or SMILES:")

if input_text:
    if all(char.isalpha() for char in input_text):
        smiles = fetch_smiles(input_text)
    else:
        smiles = input_text

    if smiles:
        with st.spinner("Generating 3D model..."):
            glb_file = generate_fake_glb(smiles)  # Replace with actual generation later
            glb_url = f"https://yourgithubusername.github.io/yourrepo/{glb_file}"  # Replace with real link

            st.subheader("📱 View in Augmented Reality")
            st.markdown(render_model_viewer(glb_url), unsafe_allow_html=True)

            st.subheader("🧪 Molecular Information")
            atom_info = get_element_info(smiles)
            for symbol, data in atom_info.items():
                st.markdown(f"**{data['name']} ({symbol})**  ")
                st.markdown(f"Atomic Number: {data['number']}  ")
                st.markdown(f"Atomic Mass: {data['mass']}  ")
                st.markdown(f"Electron Configuration: `{data['config']}`  ")
                st.markdown("---")
    else:
        st.error("Could not find SMILES for that name.")

