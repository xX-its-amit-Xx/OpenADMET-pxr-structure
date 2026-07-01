"""Replace RDKit-mangled SMILES with canonical SMILES from RCSB CCD.

The PDB→SMILES round-trip via RDKit is lossy (drops bond orders, makes phosphorus
PH instead of P, duplicates multi-copy occupancy, etc.). For Boltz-2 to predict
the correct molecule, we need the canonical SMILES from RCSB Chemical Component
Dictionary.

Reads holos.csv, queries RCSB for each lig_code, writes holos_fixed.csv + updates
boltz_inputs/*.yaml.
"""
from __future__ import annotations
import json
import time
import urllib.request
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
OUT_BASE = REPO / "data" / "processed" / "validation_set"
IN_CSV = OUT_BASE / "holos.csv"
OUT_CSV = OUT_BASE / "holos_fixed.csv"
YAML_DIR = OUT_BASE / "boltz_inputs"
PXR_FASTA = REPO.parent / "OpenADMET-pxr-challenge" / "tutorial" / "inputs" / "PXR_protein_sequence.fasta"

RCSB_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"


def fetch_smiles(lig_code: str):
    """Query RCSB ChemComp REST for canonical SMILES."""
    url = RCSB_URL.format(lig_code)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=10) as r:
            data = json.loads(r.read().decode())
    except Exception as e:
        print(f"  {lig_code}: fetch failed ({e})")
        return None
    descs = data.get("rcsb_chem_comp_descriptor", {})
    return descs.get("SMILES_stereo") or descs.get("SMILES")


def main():
    df = pd.read_csv(IN_CSV)
    pxr_seq = "".join(line.strip() for line in PXR_FASTA.read_text().splitlines() if not line.startswith(">"))
    print(f"PXR sequence: {len(pxr_seq)} residues")

    smiles_canonical = []
    bad = []
    for i, row in df.iterrows():
        code = row["lig_code"]
        print(f"[{i+1}/{len(df)}] {row['holo_id']} ({code})")
        smi = fetch_smiles(code)
        if not smi:
            bad.append(row["holo_id"])
            smiles_canonical.append(None)
        else:
            smiles_canonical.append(smi)
            print(f"  -> {smi}")
        time.sleep(0.3)  # politeness

    df["smiles_canonical"] = smiles_canonical
    df_clean = df.dropna(subset=["smiles_canonical"]).reset_index(drop=True)
    df_clean.to_csv(OUT_CSV, index=False)
    print(f"\n{len(df_clean)}/{len(df)} have canonical SMILES from RCSB")
    if bad:
        print(f"Dropped: {bad}")

    # Rewrite YAMLs with canonical SMILES
    n_written = 0
    for _, row in df_clean.iterrows():
        yaml_path = YAML_DIR / f"{row['holo_id']}.yaml"
        yaml_text = f"""sequences:
  - protein:
      id: A
      sequence: {pxr_seq}
  - ligand:
      id: L
      smiles: '{row['smiles_canonical']}'
"""
        yaml_path.write_text(yaml_text)
        n_written += 1
    print(f"Updated {n_written} YAML inputs with canonical SMILES")

    # Delete YAMLs for dropped holos
    for bid in bad:
        bad_yaml = YAML_DIR / f"{bid}.yaml"
        if bad_yaml.exists():
            bad_yaml.unlink()


if __name__ == "__main__":
    main()
