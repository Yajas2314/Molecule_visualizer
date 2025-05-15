import streamlit as st
import pubchempy as pcp
import py3dmol
from rdkit import Chem
from rdkit.Chem import Descriptors

# Get SMILES from molecule name using PubChem
def get_smiles_from_name(name):
    try:
        compounds = pcp.get_compounds(name, 'name')
        if compounds:
            return compounds[0].canonical_smiles
    except:
        pass
    return None

# Detect ionic vs covalent based on atomic composition (simple heuristic)
def detect_bond_type(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Unknown"
    # Check for presence of metals (common ionic elements)
    metals = {'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra', 'Al', 'Fe', 'Cu', 'Zn', 'Ag', 'Au', 'Hg'}
    atoms = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if atoms.intersection(metals):
        return "Ionic (contains metal ions)"
    # Check for charged groups (simplified)
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            return "Ionic (charged groups detected)"
    return "Covalent"

# Generate 3D view with style and show lone pairs approx. (highlight electronegative atoms)
def show_molecule(smiles, style='stick'):
    view = py3dmol.view(query='smiles:' + smiles)
    # Apply chosen style
    view.setStyle({style: {}})
    
    # Highlight electronegative atoms to simulate lone pairs
    electronegative = ['O', 'N', 'F', 'Cl', 'Br', 'I', 'S']
    for elem in electronegative:
        view.setStyle({'elem': elem}, {'sphere': {'color': 'yellow', 'radius': 0.3}})
    
    view.zoomTo()
    view.show()
    return view.js()

def main():
    st.title("⚗️ 3D Molecule Visualizer with Bond Type & Lone Pairs")

    mol_name = st.text_input("Enter Molecule Name (e.g. glucose):")
    smiles_code = st.text_input("Or Enter SMILES Code (e.g. C(CO)O):")
    style = st.selectbox("Choose Visualization Style:", ['stick', 'line', 'sphere', 'ball', 'cartoon'])

    smiles = None
    if mol_name:
        smiles = get_smiles_from_name(mol_name)
        if not smiles:
            st.error("Could not find molecule with that name.")
    elif smiles_code:
        smiles = smiles_code

    if smiles:
        bond_type = detect_bond_type(smiles)
        st.markdown(f"**Detected Bond Type:** {bond_type}")

        st.markdown("### 3D Structure:")
        view = show_molecule(smiles, style)
        st.components.v1.html(view, height=500)

if __name__ == "__main__":
    main()
