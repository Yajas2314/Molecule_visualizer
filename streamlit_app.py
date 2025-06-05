import streamlit as st
import py3Dmol
import requests
import re
import streamlit.components.v1 as components

# Element data for atomic numbers 1–30
element_data = {
    "H": {"name": "Hydrogen", "Z": 1, "mass": 1.008, "config": "1s1"},
    "He": {"name": "Helium", "Z": 2, "mass": 4.0026, "config": "1s2"},
    "Li": {"name": "Lithium", "Z": 3, "mass": 6.94, "config": "[He] 2s1"},
    "Be": {"name": "Beryllium", "Z": 4, "mass": 9.0122, "config": "[He] 2s2"},
    "B": {"name": "Boron", "Z": 5, "mass": 10.81, "config": "[He] 2s2 2p1"},
    "C": {"name": "Carbon", "Z": 6, "mass": 12.011, "config": "[He] 2s2 2p2"},
    "N": {"name": "Nitrogen", "Z": 7, "mass": 14.007, "config": "[He] 2s2 2p3"},
    "O": {"name": "Oxygen", "Z": 8, "mass": 15.999, "config": "[He] 2s2 2p4"},
    "F": {"name": "Fluorine", "Z": 9, "mass": 18.998, "config": "[He] 2s2 2p5"},
    "Ne": {"name": "Neon", "Z": 10, "mass": 20.180, "config": "[He] 2s2 2p6"},
    "Na": {"name": "Sodium", "Z": 11, "mass": 22.990, "config": "[Ne] 3s1"},
    "Mg": {"name": "Magnesium", "Z": 12, "mass": 24.305, "config": "[Ne] 3s2"},
    "Al": {"name": "Aluminum", "Z": 13, "mass": 26.982, "config": "[Ne] 3s2 3p1"},
    "Si": {"name": "Silicon", "Z": 14, "mass": 28.085, "config": "[Ne] 3s2 3p2"},
    "P": {"name": "Phosphorus", "Z": 15, "mass": 30.974, "config": "[Ne] 3s2 3p3"},
    "S": {"name": "Sulfur", "Z": 16, "mass": 32.06, "config": "[Ne] 3s2 3p4"},
    "Cl": {"name": "Chlorine", "Z": 17, "mass": 35.45, "config": "[Ne] 3s2 3p5"},
    "Ar": {"name": "Argon", "Z": 18, "mass": 39.948, "config": "[Ne] 3s2 3p6"},
    "K": {"name": "Potassium", "Z": 19, "mass": 39.098, "config": "[Ar] 4s1"},
    "Ca": {"name": "Calcium", "Z": 20, "mass": 40.078, "config": "[Ar] 4s2"},
    "Sc": {"name": "Scandium", "Z": 21, "mass": 44.956, "config": "[Ar] 3d1 4s2"},
    "Ti": {"name": "Titanium", "Z": 22, "mass": 47.867, "config": "[Ar] 3d2 4s2"},
    "V": {"name": "Vanadium", "Z": 23, "mass": 50.942, "config": "[Ar] 3d3 4s2"},
    "Cr": {"name": "Chromium", "Z": 24, "mass": 51.996, "config": "[Ar] 3d5 4s1"},
    "Mn": {"name": "Manganese", "Z": 25, "mass": 54.938, "config": "[Ar] 3d5 4s2"},
    "Fe": {"name": "Iron", "Z": 26, "mass": 55.845, "config": "[Ar] 3d6 4s2"},
    "Co": {"name": "Cobalt", "Z": 27, "mass": 58.933, "config": "[Ar] 3d7 4s2"},
    "Ni": {"name": "Nickel", "Z": 28, "mass": 58.693, "config": "[Ar] 3d8 4s2"},
    "Cu": {"name": "Copper", "Z": 29, "mass": 63.546, "config": "[Ar] 3d10 4s1"},
    "Zn": {"name": "Zinc", "Z": 30, "mass": 65.38, "config": "[Ar] 3d10 4s2"},
}

# Get SMILES using PubChem API
def get_smiles_from_name(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
    res = requests.get(url)
    if res.status_code == 200:
        return res.json()['PropertyTable']['Properties'][0]['IsomericSMILES']
    return None

# Display molecule using py3Dmol
def display_molecule(smiles):
    mol = py3Dmol.view(width=500, height=400)
    mol.addModel(smiles, 'smi')
    mol.setStyle({'stick': {}})
    mol.setBackgroundColor('white')
    mol.zoomTo()
    return mol

# Get atom symbols from SMILES
def get_atoms(smiles):
    return re.findall(r'[A-Z][a-z]?', smiles)

# Streamlit main app
def main():
    st.title("🧪 Molecule Visualizer (3D Viewer)")

    query = st.text_input("Enter molecule name or SMILES:", placeholder="e.g. water or CCO")
    if query:
        # Determine SMILES
        if all(c.isalnum() or c in "-=#()[]@+\\/." for c in query):
            # Likely a SMILES
            smiles = query.strip()
            st.info(f"Detected SMILES: `{smiles}`")
        else:
            # Convert name to SMILES
            smiles = get_smiles_from_name(query.strip())
            if not smiles:
                st.error("Molecule not found on PubChem.")
                return
            st.success(f"Name converted to SMILES: `{smiles}`")

        # Show molecule
        mol = display_molecule(smiles)
        components.html(mol._make_html(), height=400)

        # Element details
        atoms = set(get_atoms(smiles))
        st.subheader("🔬 Element Information")
        for atom in atoms:
            if atom in element_data:
                info = element_data[atom]
                st.markdown(f"""
                **{info['name']} ({atom})**  
                - Atomic Number: {info['Z']}  
                - Mass: {info['mass']} u  
                - Electron Configuration: `{info['config']}`
                """)
            else:
                st.markdown(f"ℹ️ No data found for atom: `{atom}`")

if __name__ == "__main__":
    main()
