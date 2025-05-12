import streamlit as st
import py3Dmol

# Set up Streamlit page
st.set_page_config(page_title="Molecule Visualizer with AR/VR", layout="centered")
st.title("🔬 Molecule Visualizer with 3D + Virtual Reality")

# Input mode
input_type = st.radio("Choose your input method:", ["Molecule Name", "SMILES Code"])

user_input = st.text_input(f"Enter the {input_type.lower()}:")

# Dictionary of molecule name to SMILES
molecule_dict = {
    "water": "O",
    "ethanol": "CCO",
    "benzene": "c1ccccc1",
    "methane": "C",
    "acetone": "CC(=O)C",
    "glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
}

# On button click
if st.button("🔍 Visualize Molecule"):
    if not user_input:
        st.error("❗ Please enter a molecule name or SMILES code.")
    else:
        if input_type == "Molecule Name":
            smiles = molecule_dict.get(user_input.lower())
            if not smiles:
                st.error("❌ Molecule name not found in database.")
                st.info("Try using a SMILES code instead.")
                st.stop()
        else:
            smiles = user_input.strip()

        # Show 3D molecule
        st.subheader("🧪 Interactive 3D Molecule Viewer")
        viewer = py3Dmol.view(width=500, height=400)
        viewer.addModel(smiles, "smi")
        viewer.setStyle({"stick": {}})
        viewer.zoomTo()
        st.components.v1.html(viewer._make_html(), height=400)

        # VR/AR section
        st.subheader("🌐 Virtual Reality Molecule View (AR Supported)")

        st.markdown("""
        <model-viewer 
            src="https://modelviewer.dev/shared-assets/models/Astronaut.glb" 
            alt="3D molecule model"
            ar 
            auto-rotate 
            camera-controls 
            style="width: 100%; height: 500px;">
        </model-viewer>

        <script type="module" src="https://unpkg.com/@google/model-viewer/dist/model-viewer.min.js"></script>
        """, unsafe_allow_html=True)

        st.caption("ℹ️ You can upload your own .glb file or replace the model URL above to show a real molecule.")
