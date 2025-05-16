import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import pyvista as pv
import tempfile
import os

# Function to generate .glb from SMILES
def generate_glb_from_smiles(smiles: str, output_file: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()

    plotter = pv.Plotter(off_screen=True)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        sphere = pv.Sphere(radius=0.3, center=(pos.x, pos.y, pos.z))

        element = atom.GetSymbol()
        color_map = {
            'H': 'white', 'C': 'black', 'O': 'red',
            'N': 'blue', 'S': 'yellow', 'Cl': 'green',
            'F': 'lime', 'Br': 'darkred', 'I': 'purple'
        }
        color = color_map.get(element, 'gray')
        plotter.add_mesh(sphere, color=color)

    # Export combined mesh
    if not plotter.mesh:
        return False
    
    plotter.mesh.save(output_file)
    return True

# Streamlit interface
st.set_page_config(page_title="Molecule to AR (.glb)", layout="centered")
st.title("🧪 Molecule to WebAR (.glb Generator)")

smiles = st.text_input("🔬 Enter SMILES string (e.g., H2O → O, CH4 → C):")

if st.button("📦 Generate .glb File"):
    if not smiles:
        st.warning("Please enter a SMILES string.")
    else:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".glb") as tmp:
            output_path = tmp.name
            success = generate_glb_from_smiles(smiles, output_path)

            if success and os.path.exists(output_path):
                with open(output_path, "rb") as f:
                    st.download_button("⬇️ Download .glb", f, file_name="molecule.glb")
                st.success("✅ .glb created! You can try it on [modelviewer.dev](https://modelviewer.dev/editor/)")
            else:
                st.error("❌ Failed to generate 3D model. Try a different SMILES.")
