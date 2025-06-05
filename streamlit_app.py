import streamlit as st
import requests

# ---------- Utility: Convert molecule name to SMILES ----------
def get_smiles_from_name(name: str) -> str | None:
    """Get SMILES string from molecule name using PubChem PUG REST API."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        smiles = data['PropertyTable']['Properties'][0]['IsomericSMILES']
        return smiles
    else:
        return None

# ---------- Utility: Render AR Viewer HTML ----------
def render_model_viewer(glb_url: str) -> str:
    """Return HTML string with <model-viewer> for AR display."""
    return f'''
    <script type="module" src="https://unpkg.com/@google/model-viewer/dist/model-viewer.min.js"></script>
    <model-viewer src="{glb_url}" alt="3D molecule model"
        ar ar-modes="webxr scene-viewer quick-look"
        auto-rotate camera-controls
        style="width: 100%; height: 500px;">
    </model-viewer>
    '''

# ---------- Main Streamlit app ----------
def main():
    st.set_page_config(page_title="AR Molecule Visualizer", layout="centered")
    st.title("🔬 AR Molecule Visualizer")

    # Input from user: SMILES or molecule name
    user_input = st.text_input("Enter molecule name or SMILES code:", placeholder="e.g. water or CCO")

    if user_input:
        # Simple heuristic: if input looks like SMILES (contains only SMILES chars)
        smiles_chars = set("CNHOPS1234567890-=#@()/\\[]")
        is_smiles = all(c in smiles_chars for c in user_input.replace(' ', ''))

        if is_smiles:
            smiles = user_input.strip()
            st.info(f"Input detected as SMILES: `{smiles}`")
        else:
            smiles = get_smiles_from_name(user_input.strip())
            if smiles is None:
                st.error("❌ Molecule name not found in PubChem database.")
                return
            st.info(f"Converted molecule name to SMILES: `{smiles}`")

        # --- Map some common SMILES to .glb files hosted on GitHub Pages ---
        # You need to upload .glb files yourself and update this dictionary accordingly
        smiles_to_glb = {
            "O": "water.glb",
            "C1=CC=CC=C1": "benzene.glb",
            "CCO": "ethanol.glb",
            # Add more mappings as you upload more .glb files
        }

        glb_filename = smiles_to_glb.get(smiles)
        if glb_filename:
            glb_url = f"https://yajas2314.github.io/molecule_visualizer/{glb_filename}"

            st.success("3D AR model available! You can view and interact below:")
            # Show the AR viewer embedded
            st.markdown(render_model_viewer(glb_url), unsafe_allow_html=True)

            # Also provide a manual "View in AR" button (optional)
            st.markdown(f'''
            <a href="{glb_url}" target="_blank" 
                style="display: inline-block; margin-top: 15px; padding: 10px 20px; background-color: #4CAF50; color: white; text-decoration: none; border-radius: 5px;">
                Open AR Model in New Tab
            </a>
            ''', unsafe_allow_html=True)

        else:
            st.warning("⚠️ 3D AR model for this molecule is not yet available.")
            st.info("You can visualize it using other molecular viewers or request the model.")

if __name__ == "__main__":
    main()
