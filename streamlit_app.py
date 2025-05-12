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
def ar_viewer_html():
    # AR Scene with A-Frame and AR.js integration
    ar_scene_html = """
    <html>
    <head>
        <script src="https://aframe.io/releases/0.9.2/aframe.min.js"></script>
        <script src="https://cdn.jsdelivr.net/gh/jeromeetienne/AR.js/aframe/build/aframe-ar.js"></script>
    </head>
    <body style="margin: 0; overflow: hidden;">
        <a-scene embedded arjs>
            <a-marker preset="hiro">
                <!-- Placeholder 3D model (box) for AR -->
                <a-entity position="0 0 0" rotation="0 0 0" scale="0.1 0.1 0.1" geometry="primitive: box; width: 1; height: 1; depth: 1;" material="color: blue;">
                </a-entity>
            </a-marker>
            <a-entity camera></a-entity>
        </a-scene>
    </body>
    </html>
    """
    return ar_scene_html

# Function to generate 3Dmol visualization of the molecule
def generate_3Dmol_visualization(smiles):
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(smiles, 'mol')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    viewer.render()

    return viewer._to_html()

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
        # Display the 3D viewer
        st.write("3D Model of the Molecule:")
        molecule_html = generate_3Dmol_visualization(smiles)
        st.components.v1.html(molecule_html, height=400)

        # Display AR Viewer with placeholder model
        st.write("AR Viewer (Placeholder Model):")
        ar_html = ar_viewer_html()
        st.components.v1.html(ar_html, height=600)
    else:
        st.warning("Could not fetch the molecule's SMILES string. Please try a different name.")
