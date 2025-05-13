# Regenerating the corrected app with properly formatted multiline HTML string

import zipfile
from pathlib import Path

# Define paths
base_dir = Path("/mnt/data/ar_molecule_streamlit_fixed")
base_dir.mkdir(exist_ok=True)
script_path = base_dir / "streamlit_app.py"
readme_path = base_dir / "README.md"

# Corrected Streamlit code with proper HTML embedding
corrected_code = '''
import streamlit as st
import requests
import py3Dmol
import tempfile
import subprocess

smiles_input = st.text_input("Enter SMILES code (optional):")
molecule_name = st.text_input("Enter molecule name:")

if smiles :
     smiles = smiles_input.strip()
elif:
     smiles = get_smiles_from_name(molecule_name)
else:
     smiles = None

st.set_page_config(layout="wide")
st.title("AR Chemistry Molecule Visualizer 🌐🧪")

def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
    res = requests.get(url)
    if res.status_code == 200:
        data = res.json()
        return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    return None

def smiles_to_sdf(smiles):
    with tempfile.NamedTemporaryFile(delete=False, suffix=".smi") as smi_file:
        smi_file.write(smiles.encode())
        smi_path = smi_file.name
    sdf_path = smi_path.replace(".smi", ".sdf")
    subprocess.run(["obabel", smi_path, "-O", sdf_path, "--gen3d"], check=True)
    return sdf_path

def sdf_to_glb(sdf_path):
    glb_path = sdf_path.replace(".sdf", ".glb")
    # Simulate the glb path (you can replace with actual conversion logic)
    return glb_path

def show_3d_molecule(smiles):
    mol = py3Dmol.view(width=400, height=400)
    mol.addModel(smiles, "smi")
    mol.setStyle({"stick": {}})
    mol.zoomTo()
    return mol

molecule = st.text_input("Enter Molecule Name (e.g., Glucose, Benzene):")

if molecule:
    with st.spinner("Generating molecule..."):
        smiles = get_smiles_from_name(molecule)
        if not smiles:
            st.error("Could not find SMILES for the molecule.")
        else:
            st.success(f"SMILES: {smiles}")
            st.subheader("3D Structure Preview")
            mol_view = show_3d_molecule(smiles)
            mol_view.show()
            try:
                sdf_path = smiles_to_sdf(smiles)
                glb_path = sdf_to_glb(sdf_path)
                st.subheader("Augmented Reality Viewer 🌐")
                ar_html = f"""
                <model-viewer 
                    src="{glb_path}" 
                    ar 
                    ar-modes="webxr scene-viewer quick-look" 
                    environment-image="neutral" 
                    auto-rotate 
                    camera-controls 
                    style="width: 100%; height: 500px;">
                </model-viewer>
                """
                st.components.v1.html(ar_html, height=520)
            except Exception as e:
                st.error(f"Error in conversion: {e}")
'''

# README content
readme_text = '''
# AR Molecule Visualizer

This app lets you:
- Input molecule name
- View SMILES and 3D structure
- Generate .glb file (simulated) for AR display using model-viewer

## Requirements
- Python 3.8+
- streamlit, requests, py3Dmol
- Open Babel (for SMILES to SDF conversion): install using conda install -c conda-forge openbabel

## Run
```bash
streamlit run streamlit_app.py
