import streamlit as st
import py3Dmol
import requests

# Function to fetch SMILES from PubChem based on molecule name
def get_smiles_from_pubchem(molecule_name):
    # PubChem API to search for the molecule by name
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{molecule_name}/property/SMILES/TXT"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        return response.text.strip()  # Return the SMILES string
    except requests.exceptions.RequestException as e:
        st.error(f"Error fetching SMILES for {molecule_name}: {e}")
        return ""

# Function to generate AR scene HTML with A-Frame and AR.js
def ar_viewer_html(smiles):
    # Generate the molecule using py3Dmol
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(smiles, 'mol')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    viewer.render()

    # Get the HTML representation of the 3Dmol model
    viewer_str = viewer._to_html()

    # AR Scene with A-Frame and AR.js integration
    ar_scene_html = f"""
    <html>
    <head>
        <script src="https://aframe.io/releases/0.9.2/aframe.min.js"></script>
        <script src="https://cdn.jsdelivr.net/gh/jeromeetienne/AR.js/aframe/build/aframe-ar.js"></script>
    </head>
    <body style="margin: 0; overflow: hidden;">
        <a-scene embedded arjs>
            <a-marker preset="hiro">
                <a-entity position="0 0 0" rotation="0 0 0" scale="0.1 0.1 0.1" geometry="primitive: box; width: 1; height: 1; depth: 1;" material="color: blue;">
                </a-entity>
            </a-marker>
            <a-entity camera></a-entity>
        </a-scene>
    </body>
    </html>
    """
    return ar_scene_html

# Streamlit App UI
st.title("Molecular AR Viewer")
st.write("""
         This app allows you to view molecular structures in Augmented Reality (AR). 
         Enter a molecule name (like water, benzene, methane) to display its structure in AR.
         """)

# Input: Molecule Name
molecule_name = st.text_input("Enter Molecule Name", "")

if molecule_name:
    # Fetch the SMILES string from PubChem API
    smiles = get_smiles_from_pubchem(molecule_name)

    if smiles:
        # Display the AR Viewer with the corresponding molecule
        ar_html = ar_viewer_html(smiles)
        st.components.v1.html(ar_html, height=600)
    else:
        st.warning("Could not fetch the molecule's SMILES string. Please try a different name.")
