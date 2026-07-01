"""Build validation set from cached PXR holo PDBs.

For each holo PDB:
  1. Extract the bound ligand: SMILES (canonical) + crystal HETATM coords
  2. Skip if no drug-like ligand (filter by MW, exclude solvent codes)
  3. Save to:
     - data/processed/validation_set/holos.csv: holo_id,smiles,lig_code,mw,n_atoms
     - data/processed/validation_set/crystal_ligands/{holo_id}.pdb: just the LIG block
     - data/processed/validation_set/boltz_inputs/{holo_id}.yaml: Boltz-2 input

Boltz-2 input format (yaml):
    sequences:
      - protein:
          id: A
          sequence: <PXR FASTA>
      - ligand:
          id: L
          smiles: <ligand SMILES>
"""
from __future__ import annotations
import shutil
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

REPO = Path(__file__).resolve().parent.parent
HOLO_DIRS = [REPO / "data" / "external" / "pxr_holo",
             REPO / "data" / "external" / "pxr_crystals"]
PXR_FASTA = REPO.parent / "OpenADMET-pxr-challenge" / "tutorial" / "inputs" / "PXR_protein_sequence.fasta"
OUT_BASE = REPO / "data" / "processed" / "validation_set"
OUT_CSV = OUT_BASE / "holos.csv"
OUT_LIG_DIR = OUT_BASE / "crystal_ligands"
OUT_YAML_DIR = OUT_BASE / "boltz_inputs"

SKIP_RES = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "MES", "TRS", "DMS",
            "NAG", "PEG", "CIT", "FMT", "ACT", "MPD", "PGE", "DOD",
            "CA", "MG", "ZN", "CL", "NA", "K", "NO3", "NHE", "IMD", "IPA",
            "BME", "DTT", "TAR", "EPE", "PEP", "NH4", "BR", "IOD", "GTP",
            "ANP", "HEM", "NDP"}
MW_MIN = 150


def extract_ligand_records(pdb_path):
    """Return list of (lig_code, [HETATM lines], [CONECT lines for those serials])."""
    text = pdb_path.read_text(errors="replace")
    lines = text.splitlines()

    # Map lig_code → atom serials, hetatm lines
    code_to_serials = {}
    code_to_lines = {}
    for line in lines:
        if line.startswith("HETATM"):
            res = line[17:20].strip()
            if res in SKIP_RES or len(res) < 2:
                continue
            try:
                serial = int(line[6:11])
            except ValueError:
                continue
            code_to_serials.setdefault(res, set()).add(serial)
            code_to_lines.setdefault(res, []).append(line)

    out = []
    for code, serials in code_to_serials.items():
        het_lines = code_to_lines[code]
        if len(het_lines) < 8:  # minimum atoms for "real" ligand
            continue
        conect_lines = []
        for line in lines:
            if line.startswith("CONECT"):
                parts = line.split()
                if any(p.isdigit() and int(p) in serials for p in parts[1:]):
                    conect_lines.append(line)
        out.append((code, het_lines, conect_lines))
    return out


def lig_lines_to_mol(het_lines, conect_lines):
    mini = "REMARK 1\n" + "\n".join(het_lines + conect_lines) + "\nEND\n"
    try:
        mol = Chem.MolFromPDBBlock(mini, sanitize=True, removeHs=True)
        if mol and mol.GetNumAtoms() >= 8:
            return mol
    except Exception:
        return None
    return None


def main():
    if OUT_BASE.exists():
        shutil.rmtree(OUT_BASE)
    OUT_LIG_DIR.mkdir(parents=True)
    OUT_YAML_DIR.mkdir(parents=True)

    # Load PXR FASTA
    if not PXR_FASTA.exists():
        print(f"ERROR: PXR FASTA missing at {PXR_FASTA}")
        return
    fasta_text = PXR_FASTA.read_text()
    pxr_seq = "".join(line.strip() for line in fasta_text.splitlines() if not line.startswith(">"))
    print(f"PXR sequence: {len(pxr_seq)} residues")

    records = []
    seen_smiles = set()

    for holo_dir in HOLO_DIRS:
        if not holo_dir.exists():
            continue
        for pdb in sorted(holo_dir.glob("*.pdb")):
            holo_id = pdb.stem
            ligs = extract_ligand_records(pdb)
            if not ligs:
                continue
            # Use the first valid ligand
            for code, het_lines, conect_lines in ligs:
                mol = lig_lines_to_mol(het_lines, conect_lines)
                if mol is None:
                    continue
                try:
                    smi = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
                    mw = Descriptors.MolWt(mol)
                    n_atoms = mol.GetNumAtoms()
                except Exception:
                    continue
                if mw < MW_MIN:
                    continue
                if smi in seen_smiles:
                    print(f"  Skip {holo_id} ({code}) — duplicate SMILES already in set")
                    break
                seen_smiles.add(smi)
                # Write the LIG-only PDB
                lig_pdb_path = OUT_LIG_DIR / f"{holo_id}.pdb"
                lig_pdb_path.write_text("REMARK 1 crystal ligand from " + holo_id + "\n" +
                                         "\n".join(het_lines + conect_lines) + "\nEND\n")
                # Write Boltz-2 input YAML
                yaml_path = OUT_YAML_DIR / f"{holo_id}.yaml"
                yaml_text = f"""sequences:
  - protein:
      id: A
      sequence: {pxr_seq}
  - ligand:
      id: L
      smiles: '{smi}'
"""
                yaml_path.write_text(yaml_text)
                records.append({
                    "holo_id": holo_id,
                    "smiles": smi,
                    "lig_code": code,
                    "mw": round(mw, 2),
                    "n_atoms": n_atoms,
                })
                break  # one ligand per holo
    df = pd.DataFrame(records)
    df.to_csv(OUT_CSV, index=False)
    print(f"\nValidation set: {len(df)} holos with valid ligands")
    print(f"  CSV: {OUT_CSV}")
    print(f"  Crystal ligands: {OUT_LIG_DIR}")
    print(f"  Boltz inputs: {OUT_YAML_DIR}")
    print(f"\nMW range: {df['mw'].min():.0f} – {df['mw'].max():.0f}")
    print(f"N_atoms range: {df['n_atoms'].min()} – {df['n_atoms'].max()}")


if __name__ == "__main__":
    main()
