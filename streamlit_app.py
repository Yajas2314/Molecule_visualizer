import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import requests
import periodictable

st.set_page_config(page_title="Advanced Molecule Visualizer", layout="wide")
st.title("🧪 Molecule Visualizer with Lone Pairs & Configurations")

# Input
user_input = st.text_input("🔍 Enter molecule name or SMILES")

def fetch_smiles_from_name(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        res = requests.get(url)
        data = res.json()
        return data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
    except:
        return None

def get_molecule_details(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/JSON"
        res = requests.get(url)
        data = res.json()
        info = data["PC_Compounds"][0]
        props = info.get("props", [])
        weight = next((p["value"]["fval"] for p in props if p.get("urn", {}).get("label") == "Molecular Weight"), "N/A")
        formula = next((p["value"]["sval"] for p in props if p.get("urn", {}).get("label") == "Molecular Formula"), "N/A")
        atoms = info.get("atoms", {}).get("element", [])
        return {"Formula": formula, "Weight": weight, "Atoms": atoms}
    except:
        return None

def get_e_config(symbol):
    try:
        element = getattr(periodictable, symbol.lower())
        config = ""
        electrons = element.number
        orbitals = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s']
        capacity = {'s': 2, 'p': 6, 'd': 10}
        for orb in orbitals:
            if electrons <= 0:
                break
            subshell = orb[-1]
            max_e = capacity[subshell]
            fill = min(electrons, max_e)
            config += f"{orb}{fill} "
            electrons -= fill
        return config.strip()
    except:
        return "Not found"

if user_input:
    if all(c.isalnum() or c in "=#@[]()\\/.-+" for c in user_input):
        smiles = user_input
        molname = None
    else:
        molname = user_input
        smiles = fetch_smiles_from_name(user_input)

    if smiles:
        st.subheader("🔬 3D Molecular Structure")
        view = py3Dmol.view(width=600, height=450)
        view.addModel(smiles, "smi")
        view.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
        view.setBackgroundColor("white")
        view.zoomTo()
        components.html(view.render().data, height=470)

        st.subheader("📘 Molecule Details")
        details = get_molecule_details(molname or smiles)
        if details:
            st.write(f"**Formula:** {details['Formula']}")
            st.write(f"**Molecular Weight:** {details['Weight']} g/mol")
        else:
            st.info("No compound data available.")

        if details and details['Atoms']:
            st.subheader("🔬 Electronic Configuration of Elements")
            unique_atoms = set(details["Atoms"])
            for atomic_num in unique_atoms:
                try:
                    symbol = periodictable.elements[atomic_num].symbol
                    config = get_e_config(symbol)
                    st.markdown(f"**{symbol} (Z={atomic_num})** → {config}")
                except:
                    continue

        st.subheader("💡 Note on Lone Pairs")
        st.markdown("""
        - Lone pairs are not directly encoded in SMILES.
        - The structure uses **valence rules**; lone pairs are **approximated visually** via shape.
        - You can infer lone pairs based on the **valence shell electron pair repulsion (VSEPR)** model.
        """)
    else:
        st.error("❌ Invalid molecule name or SMILES.")
