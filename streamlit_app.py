import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import requests
import periodictable
import hashlib

st.set_page_config(page_title="AR Molecule Visualizer", layout="wide")
st.title("🧪 AR/VR Molecule Visualizer")

user_input = st.text_input("🔍 Enter molecule name or SMILES")

def fetch_smiles(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        r = requests.get(url)
        return r.json()['PropertyTable']['Properties'][0]['IsomericSMILES']
    except:
        return None

def fetch_info(name):
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/JSON"
        r = requests.get(url)
        data = r.json()
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
        el = getattr(periodictable, symbol.lower())
        z = el.number
        orbitals = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s']
        fill = {'s': 2, 'p': 6, 'd': 10}
        config = ""
        for orb in orbitals:
            if z <= 0: break
            kind = orb[-1]
            e = min(fill[kind], z)
            config += f"{orb}{e} "
            z -= e
        return config
    except:
        return "Unknown"

def get_ar_url(smiles):
    # Use MD5 hash of SMILES as filename
    code = hashlib.md5(smiles.encode()).hexdigest()
    return f"https://raw.githubusercontent.com/Yajas2314/Molecule_visualizer/main/models/{code}.glb"

if user_input:
    # Convert to SMILES if needed
    smiles = user_input if all(c.isalnum() or c in "=#@[]()\\/.-+" for c in user_input) else fetch_smiles(user_input)

    if smiles:
        st.subheader("🧬 3D Molecular View")
        viewer = py3Dmol.view(width=600, height=450)
        viewer.addModel(smiles, "smi")
        viewer.setStyle({'stick': {}, 'sphere': {'scale': 0.3}})
        viewer.setBackgroundColor("white")
        viewer.zoomTo()

        # ✅ FIXED HERE:
        components.html(viewer._make_html(), height=470)

        # AR model viewer
        st.subheader("👓 View in Augmented Reality")
        ar_url = get_ar_url(smiles)
        st.markdown(f"""
        <model-viewer src="{ar_url}" alt="Molecule" ar ar-modes="scene-viewer webxr quick-look" auto-rotate camera-controls 
        style="width: 100%; height: 400px;"></model-viewer>

        <script type="module" src="https://unpkg.com/@google/model-viewer/dist/model-viewer.min.js"></script>
        """, unsafe_allow_html=True)

        # Details
        st.subheader("📘 Molecular Information")
        info = fetch_info(user_input)
        if info:
            st.write(f"**Formula**: {info['Formula']}")
            st.write(f"**Molecular Weight**: {info['Weight']} g/mol")
            st.subheader("🧠 Electronic Configurations")
            for z in set(info["Atoms"]):
                try:
                    symbol = periodictable.elements[z].symbol
                    st.markdown(f"**{symbol}**: {get_e_config(symbol)}")
                except:
                    continue
        else:
            st.info("No extra info available.")
    else:
        st.error("Invalid molecule name or SMILES format.")
