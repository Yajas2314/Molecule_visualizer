import streamlit as st
import requests
from PIL import Image
from io import BytesIO

st.set_page_config(page_title="Molecule Visualizer", layout="centered")
st.title("🧪 Molecule Structure Viewer")
st.write("Enter a valid SMILES code (e.g., CCO for ethanol):")

smiles = st.text_input("SMILES Code", "CCO")

if smiles:
    try:
        # Use NCI’s Cactus API to get the structure image
        url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/image"
        response = requests.get(url)

        if response.status_code == 200:
            img = Image.open(BytesIO(response.content))
            st.image(img, caption=f"Structure of {smiles}")
        else:
            st.error("❌ Invalid SMILES or unable to fetch image.")
    except Exception as e:
        st.error(f"⚠️ Error: {e}")
