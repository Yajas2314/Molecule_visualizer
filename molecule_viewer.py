import streamlit as st
import py3Dmol

# Streamlit App Title
st.set_page_config(page_title="3D Molecule Viewer", layout="wide")
st.title("🔬 3D Molecule Viewer")
st.write("Enter a SMILES code to visualize the 3D structure.")

# Input from user
smiles_input = st.text_input("Enter SMILES code (e.g. CCO for ethanol):", "CCO")

# Only show viewer if SMILES is entered
if smiles_input:
    try:
        # Create 3Dmol.js viewer
        viewer = py3Dmol.view(width=600, height=400)
        viewer.addModel(smiles_input, "smi")  # Interpret as SMILES
        viewer.setStyle({"stick": {}})
        viewer.zoomTo()

        # Render in Streamlit
        viewer_html = viewer._make_html()
        st.components.v1.html(viewer_html, height=400)
    except Exception as e:
        st.error(f"❌ Error rendering molecule: {e}")
