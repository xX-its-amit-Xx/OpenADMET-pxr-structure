"""Post-processing helper for nb16 results.

When nb16 finishes on Kaggle, run:
  python scripts/kaggle_push.py --nb 16 --pull
to download outputs to submissions/kaggle_nb16/.

Then run THIS script to:
  1) Extract HETATM LIG records from each of the 18 new crystal PDBs into
     data/processed/validation_set/crystal_ligands/<PDBID>.pdb
  2) Unzip the Kaggle output zip into a predictions dir laid out like nb06
     (one subdir per holo containing the 20 *.pdb predictions)
  3) Tell the caller to run score_validation.py to extend scores.csv.

Reads:
  data/processed/chembl_pxr_ranker/new_validation_pdbs.csv
  submissions/kaggle_nb16/*.zip                 (Boltz outputs)

Writes:
  data/processed/validation_set/crystal_ligands/<PDBID>.pdb   (ligand-only)
  data/processed/validation_set/preds_nb16/<HOLO>/*.pdb
"""
from __future__ import annotations
import sys, zipfile, shutil
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
NEW_CSV = REPO / "data" / "processed" / "chembl_pxr_ranker" / "new_validation_pdbs.csv"
VAL_DIR = REPO / "data" / "processed" / "validation_set"
CRY_DIR = VAL_DIR / "crystal_ligands"
PRED_DIR = VAL_DIR / "preds_nb16"
KAGGLE_OUT = REPO / "submissions" / "kaggle_nb16"


def extract_lig(pdb_path: Path, lig_code: str, out_path: Path) -> bool:
    """Write only HETATM lines whose residue name matches lig_code (or first 3 chars)."""
    keep = []
    code = lig_code.strip().upper()[:3]
    for line in pdb_path.read_text(errors="replace").splitlines():
        if line.startswith("HETATM"):
            res = line[17:20].strip().upper()
            if res == code:
                keep.append(line)
    if not keep:
        return False
    out_path.write_text("\n".join(keep) + "\nEND\n")
    return True


def main():
    df = pd.read_csv(NEW_CSV)
    CRY_DIR.mkdir(parents=True, exist_ok=True)
    PRED_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/3] Extracting crystal ligands for 18 new PDBs...")
    ok = 0
    for r in df.itertuples():
        src = Path(r.pdb_path)
        out = CRY_DIR / f"{r.pdb_id}.pdb"
        if out.exists():
            ok += 1
            continue
        if not src.exists():
            print(f"  MISSING source PDB: {r.pdb_id} -> {src}")
            continue
        if extract_lig(src, r.chem_comp_id, out):
            ok += 1
        else:
            print(f"  no LIG extracted for {r.pdb_id} (chem_comp={r.chem_comp_id})")
    print(f"  -> {ok}/{len(df)} crystal ligands ready in {CRY_DIR}")

    print("[2/3] Unzipping Kaggle nb16 output...")
    zips = list(KAGGLE_OUT.glob("*.zip"))
    if not zips:
        print(f"  No zips in {KAGGLE_OUT}; run kaggle_push.py --nb 16 --pull first.")
        sys.exit(1)
    for z in zips:
        print(f"  unzip {z.name}")
        with zipfile.ZipFile(z) as zf:
            zf.extractall(PRED_DIR)
    # Boltz default output layout: boltz_results_<yaml_stem>/predictions/<stem>/<stem>_model_*.pdb
    # Re-organise into PRED_DIR/<HOLO>/*.pdb so score_validation.py can pick it up.
    print("[3/3] Re-organising predictions into per-holo dirs...")
    moved = 0
    for r in df.itertuples():
        holo = r.pdb_id
        target = PRED_DIR / holo
        target.mkdir(exist_ok=True)
        # Look for any boltz_results dir matching this holo's YAML stem
        for cand in PRED_DIR.rglob("*.pdb"):
            if holo.lower() in cand.as_posix().lower() and cand.parent != target:
                shutil.copy(cand, target / cand.name)
                moved += 1
    print(f"  -> moved {moved} prediction PDBs into per-holo dirs")
    print()
    print("Done. Next:")
    print("  python scripts/score_validation.py data/processed/validation_set/preds_nb16")
    print("  python scripts/aggregate_validation.py")


if __name__ == "__main__":
    main()
