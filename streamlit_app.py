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
        return None

# Function to generate AR scene HTML with A-Frame and AR.js
def ar_viewer_html_with_model(model_url):
    ar_scene_html = f"""
    <html>
    <head>
        <script src="https://aframe.io/releases/1.2.0/aframe.min.js"></script>
        <script src="https://cdn.jsdelivr.net/gh/AR-js-org/AR.js@3.3.2/aframe/build/aframe-ar.min.js"></script>
    </head>
    <body style="margin: 0; overflow: hidden;">
        <a-scene embedded arjs>
            <a-marker preset="hiro">
                <a-entity gltf-model="{model_url}" scale="0.5 0.5 0.5" rotation="0 0 0" position="0 0 0"></a-entity>
            </a-marker>
            <a-entity camera></a-entity>
        </a-scene>
    </body>
    </html>
    """
    return ar_scene_html

# Function to generate 3Dmol visualization of the molecule
def generate_3Dmol_visualization(smiles):
    viewer = py3Dmol.view(width=800, height=400)
    viewer.addModel(smiles, 'mol')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    viewer.render()

# Get the HTML content of the viewer
viewer_html = viewer._make_html()
return viewer_html

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
        
        # Generate the 3Dmol visualization
        molecule_html = generate_3Dmol_visualization(smiles)

        # Use Streamlit to display the 3Dmol HTML
        st.components.v1.html(molecule_html, height=400)

        # Display AR Viewer with placeholder model
        st.write("AR Viewer (Placeholder Model):")
        
    def ar_viewer_html_with_model(model_url):
        return f"""
        <html>
        <head>
        <script src="https://aframe.io/releases/1.2.0/aframe.min.js"></script>
        <script src="https://cdn.jsdelivr.net/gh/AR-js-org/AR.js@3.3.2/aframe/build/aframe-ar.min.js"></script>
        </head>
        <body style="margin: 0; overflow: hidden;">
        <a-scene embedded arjs>
            <a-marker preset="hiro">
                <a-entity gltf-model="{model_url}" scale="0.5 0.5 0.5" rotation="0 0 0" position="0 0 0"></a-entity>
            </a-marker>
            <a-entity camera></a-entity>
        </a-scene>
    </body>
    </html>
    """
        
        ar_html = ar_viewer_html()

        # Use Streamlit to display AR viewer HTML
        st.components.v1.html(ar_html, height=600)
else:
          st.warning("Could not fetch the molecule's SMILES string. Please try a different name.")
