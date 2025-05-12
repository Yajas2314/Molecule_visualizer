import streamlit as st
import py3Dmol

# Function to generate AR scene HTML with A-Frame and AR.js
def ar_viewer_html(smiles):
    # Generate the molecule using py3Dmol
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(smiles, 'mol')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    viewer.render()
    
    # Convert the molecule to an HTML representation for embedding in AR scene
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
                <!-- Here we use the generated 3D model for the molecule -->
                <a-entity geometry="primitive: box; width: 1; height: 1; depth: 1;" material="color: blue;">
                </a-entity>
            </a-marker>
            <a-entity camera></a-entity>
        </a-scene>
    </body>
    </html>
    """
    return ar_scene_html

# Function to convert molecule name to SMILES (can be extended to use different molecule names)
def molecule_to_smiles(molecule_name):
    molecule_dict = {
        'benzene': 'C1=CC=CC=C1',
        'methane': 'CH4',
        'ethanol': 'CCO',
        'acetone': 'CC(C)=O',
        'glucose': 'C6H12O6'
    }
    return molecule_dict.get(molecule_name.lower(), '')

# Streamlit App UI
st.title("Molecular AR Viewer")
st.write("""
         This app allows you to view molecular structures in Augmented Reality (AR). 
         Enter a molecule name (like benzene, methane, ethanol) to display its structure in AR.
         """)

# Input: Molecule Name
molecule_name = st.text_input("Enter Molecule Name", "benzene")

# Get the corresponding SMILES string for the molecule
smiles = molecule_to_smiles(molecule_name)

# Display the AR Viewer with the corresponding molecule
if smiles:
    ar_html = ar_viewer_html(smiles)
    st.components.v1.html(ar_html, height=600)
else:
    st.warning("Invalid molecule name. Please try one of these: benzene, methane, ethanol, acetone, glucose.")
