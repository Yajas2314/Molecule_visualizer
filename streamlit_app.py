import streamlit as st
import requests
from PIL import Image
from io import BytesIO

# Function to get SMILES code from molecule name
def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        try:
            smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
            return smiles
        except KeyError:
            return None
    return None

# Streamlit UI setup
st.set_page_config(page_title="Molecule Visualizer", layout="centered")
st.title("🧪 Molecule Structure Viewer")

# Option for SMILES code or molecule name
option = st.radio("Choose input type", ["SMILES Code", "Molecule Name"])

if option == "SMILES Code":
    smiles = st.text_input("Enter SMILES Code (e.g., CCO for Ethanol)", "CCO")
    molecule_name = ""
else:
    molecule_name = st.text_input("Enter Molecule Name (e.g., Ethanol)")
    smiles = ""

if smiles or molecule_name:
    try:
        # If molecule name is provided, fetch its SMILES code
        if molecule_name:
            smiles = get_smiles_from_name(molecule_name)
            if smiles is None:
                st.error(f"❌ Could not find SMILES for {molecule_name}.")
                st.stop()
            else:
                st.write(f"Found SMILES for {molecule_name}: {smiles}")

        # Fetch molecule image from NCI's Cactus API
        url = f"https://cactus.nci.nih.gov/chemical/structure/{smiles}/image"
        response = requests.get(url)

        if response.status_code == 200:
            img = Image.open(BytesIO(response.content))
            st.image(img, caption=f"Structure of {smiles}")
        else:
            st.error("❌ Invalid SMILES or unable to fetch image.")
    except Exception as e:
        st.error(f"⚠️ Error: {e}")
