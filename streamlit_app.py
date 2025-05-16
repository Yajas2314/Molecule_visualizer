import streamlit as st
import py3Dmol
import requests
import streamlit.components.v1 as components

# Function to get SMILES from molecule name
def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/TXT"
    response = requests.get(url)
    return response.text.strip() if response.status_code == 200 else None

# Function to generate 3D viewer and return HTML
def generate_html_3d(smiles):
    mol_view = py3Dmol.view(width=600, height=400)
    mol_view.addModel(smiles, "smi")
    mol_view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
    mol_view.setBackgroundColor('white')
    mol_view.zoomTo()
    return mol_view._make_html()

# Streamlit UI
st.title("3D Molecule Viewer and HTML Export")

user_input = st.text_input("Enter molecule name or SMILES code")

if user_input:
    if all(char.isalpha() or char.isdigit() for char in user_input.replace("-", "")):
        smiles = get_smiles_from_name(user_input)
        if smiles:
            st.success(f"Found SMILES: {smiles}")
        else:
            st.error("Could not fetch SMILES for the name.")
            st.stop()
    else:
        smiles = user_input

    html_code = generate_html_3d(smiles)

    # Display the 3D molecule in Streamlit
    components.html(html_code, height=450, width=700)

    # Option to download HTML
    st.download_button("Download as HTML", data=html_code, file_name="molecule.html", mime="text/html")
