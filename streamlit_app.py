import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io

st.set_page_config(page_title="3D Molecule Viewer", layout="centered")
st.title("🔬 3D Molecule Visualizer")
st.write("Enter a valid SMILES code to see its molecular structure:")

smiles = st.text_input("SMILES Code", "CCO")  # Example: Ethanol

try:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        st.image(img, caption="Molecular Structure", use_column_width=False)
    else:
        st.error("❌ Invalid SMILES code. Please enter a correct one.")
except Exception as e:
    st.error(f"⚠️ Error: {str(e)}")

   
