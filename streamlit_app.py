import streamlit as st
import py3dmol
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem

def get_smiles_from_name(name):
    try:
        compounds = pcp.get_compounds(name, 'name')
        if compounds:
            return compounds[0].canonical_smiles
    except:
        return None

def add_hydrogens(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # Add explicit hydrogens
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)  # 3D coords
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.MolToMolBlock(mol)

def detect_ionic_or_covalent(smiles):
    ionic_elements = ['Na', 'K', 'Cl', 'Br', 'I', 'Ca', 'Mg', 'Li', 'F']
    mol = Chem.MolFromSmiles(smiles)
    atoms = [mol.GetAtomWithIdx(i) for i in range(mol.GetNumAtoms())]
    elements = set([atom.GetSymbol() for atom in atoms])
    # Simple heuristic:
    if any(e in ionic_elements for e in elements):
        return "Likely Ionic Compound"
    else:
        return "Likely Covalent Compound"

def get_lone_pairs_info(smiles):
    # Simple rough estimation by valence - bonds, not 100% accurate
    mol = Chem.MolFromSmiles(smiles)
    lone_pairs = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        valence = Chem.GetPeriodicTable().GetDefaultValence(symbol)
        bonds = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        lone_pair_count = valence - bonds
        if lone_pair_count > 0:
            lone_pairs.append(f"{symbol}{atom.GetIdx()} has approx {int(lone_pair_count)} lone pair(s)")
    return lone_pairs

def show_molecule(molblock, style='stick'):
    view = py3dmol.view(width=500, height=400)
    view.addModel(molblock, 'mol')
    view.setStyle({style: {}})
    view.zoomTo()
    html = view._make_html()
    st.components.v1.html(html, height=450)

def main():
    st.title("🔬 Advanced 3D Molecule Visualizer")
    mol_name = st.text_input("Enter Molecule Name (e.g. glucose):")
    smiles_code = st.text_input("Or Enter SMILES Code (e.g. C(CO)O):")
    style = st.selectbox("Choose Visualization Style:", ['stick', 'line', 'sphere', 'ball'])

    smiles = None
    if mol_name:
        smiles = get_smiles_from_name(mol_name)
        if not smiles:
            st.error("Molecule name not found.")
    elif smiles_code:
        smiles = smiles_code

    if smiles:
        st.write(f"**Input SMILES:** {smiles}")

        molblock = add_hydrogens(smiles)
        show_molecule(molblock, style)

        compound_type = detect_ionic_or_covalent(smiles)
        st.info(compound_type)

        lone_pairs = get_lone_pairs_info(smiles)
        if lone_pairs:
            st.write("### Estimated Lone Pairs on Atoms:")
            for lp in lone_pairs:
                st.write(f"- {lp}")
        else:
            st.write("No lone pairs detected or molecule data insufficient.")

if __name__ == "__main__":
    main()
