"""Local LDDT-PLI proxy scoring for validation set.

Given predicted complex PDBs (e.g., from Boltz-2 on the 35 holo SMILES) and the
crystal LIG PDB, compute:
  - lig_rmsd: symmetry-corrected RMSD between predicted ligand heavy atoms and
              crystal ligand heavy atoms (after protein Cα superposition)
  - lddt_pli_proxy: 1 / (1 + lig_rmsd) so lower-is-worse becomes 0-1 like real LDDT-PLI

For each validation holo, accumulates per-prediction scores so re-ranking methods
can be evaluated as "if method picks pose X, expected LDDT-PLI = ?"

Usage:
  score_validation.py PREDICTIONS_DIR
where PREDICTIONS_DIR has subdirs per holo (e.g., 1ILH/, 2CHW/) each containing
*.pdb predicted complexes. Output: data/processed/validation_set/scores.csv
"""
from __future__ import annotations
import sys, csv
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

REPO    = Path(__file__).resolve().parent.parent
VAL_DIR = REPO / "data" / "processed" / "validation_set"
CRYSTAL = VAL_DIR / "crystal_ligands"
OUT_CSV = VAL_DIR / "scores.csv"


def parse_pdb_records(text):
    """Return (atom_records, hetatm_records, conect_records).

    Each atom/hetatm record is dict(serial, name, res, chain, resnum, x, y, z, elem).
    """
    atoms, hetatms, conects = [], [], []
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "ATOM\t")):
            try:
                atoms.append({
                    "serial": int(line[6:11]),
                    "name": line[12:16].strip(),
                    "res": line[17:20].strip(),
                    "chain": line[21:22].strip(),
                    "resnum": int(line[22:26].strip() or 0),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "elem": line[76:78].strip() or line[12:14].strip(),
                    "raw": line,
                })
            except (ValueError, IndexError):
                continue
        elif line.startswith("HETATM"):
            try:
                hetatms.append({
                    "serial": int(line[6:11]),
                    "name": line[12:16].strip(),
                    "res": line[17:20].strip(),
                    "chain": line[21:22].strip(),
                    "resnum": int(line[22:26].strip() or 0),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "elem": line[76:78].strip() or line[12:14].strip(),
                    "raw": line,
                })
            except (ValueError, IndexError):
                continue
        elif line.startswith("CONECT"):
            conects.append(line)
    return atoms, hetatms, conects


def lig_from_records(hetatms, conects, lig_resname=None):
    """Build an RDKit mol from a subset of HETATM records (LIG by default)."""
    if lig_resname:
        sel = [h for h in hetatms if h["res"] == lig_resname]
    else:
        sel = hetatms
    if not sel:
        return None
    serials = {h["serial"] for h in sel}
    lines = [h["raw"] for h in sel]
    rel_conects = []
    for c in conects:
        parts = c.split()
        if any(p.isdigit() and int(p) in serials for p in parts[1:]):
            rel_conects.append(c)
    mini = "REMARK 1\n" + "\n".join(lines + rel_conects) + "\nEND\n"
    try:
        mol = Chem.MolFromPDBBlock(mini, sanitize=True, removeHs=True)
        return mol if (mol and mol.GetNumAtoms() >= 4) else None
    except Exception:
        return None


def kabsch_transform(src, dst):
    """Return (R, t) so that R @ src.T + t aligns src points to dst points."""
    src_c = src.mean(axis=0)
    dst_c = dst.mean(axis=0)
    src0 = src - src_c
    dst0 = dst - dst_c
    H = src0.T @ dst0
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = dst_c - R @ src_c
    return R, t


def superpose_proteins(pred_atoms, ref_atoms):
    """Compute (R, t) aligning predicted Cα to reference Cα.

    Ref Cα subset = whatever residues appear in both. Match by resnum (since pred
    typically starts at 1 but ref starts at the PDB residue number).
    """
    pred_ca = {a["resnum"]: np.array([a["x"], a["y"], a["z"]]) for a in pred_atoms if a["name"] == "CA"}
    ref_ca = {a["resnum"]: np.array([a["x"], a["y"], a["z"]]) for a in ref_atoms if a["name"] == "CA"}
    if not pred_ca or not ref_ca:
        return None, None, 0
    # Try direct match first
    common = sorted(set(pred_ca) & set(ref_ca))
    if len(common) >= 30:
        src = np.array([pred_ca[k] for k in common])
        dst = np.array([ref_ca[k] for k in common])
        R, t = kabsch_transform(src, dst)
        return R, t, len(common)
    # No matching resnums — try sequential pairing by sorted order (Boltz output renumbers)
    pred_keys = sorted(pred_ca.keys())
    ref_keys = sorted(ref_ca.keys())
    if not pred_keys or not ref_keys:
        return None, None, 0
    # Pair shortest to longest from start, accept some offset
    n = min(len(pred_keys), len(ref_keys))
    if n < 30:
        return None, None, 0
    src = np.array([pred_ca[k] for k in pred_keys[:n]])
    dst = np.array([ref_ca[k] for k in ref_keys[:n]])
    R, t = kabsch_transform(src, dst)
    return R, t, n


def transform_coords(coords, R, t):
    return (R @ coords.T).T + t


def sym_corrected_rmsd(pred_mol, crystal_mol):
    """RDKit's GetBestRMS does symmetry-aware RMSD on heavy atoms."""
    try:
        # Both mols must have same atoms (we use heavy-atom canonical mol)
        from rdkit.Chem import rdMolAlign
        # Use GetBestRMS — handles symmetry like phenyl ring equivalence
        return rdMolAlign.GetBestRMS(pred_mol, crystal_mol)
    except Exception:
        # Fall back to simple atom-index RMSD
        try:
            n = min(pred_mol.GetNumAtoms(), crystal_mol.GetNumAtoms())
            pc = pred_mol.GetConformer()
            cc = crystal_mol.GetConformer()
            sq = 0.0
            for i in range(n):
                pp = pc.GetAtomPosition(i)
                qq = cc.GetAtomPosition(i)
                sq += (pp.x - qq.x) ** 2 + (pp.y - qq.y) ** 2 + (pp.z - qq.z) ** 2
            return float(np.sqrt(sq / n))
        except Exception:
            return None


def score_one(pred_path, crystal_lig_mol, crystal_protein_atoms):
    """Score a single predicted complex PDB against crystal."""
    text = pred_path.read_text(errors="replace")
    atoms, hetatms, conects = parse_pdb_records(text)
    pred_lig_mol = lig_from_records(hetatms, conects, lig_resname="LIG")
    if pred_lig_mol is None:
        return None
    # Superpose protein
    R, t, n_ca = superpose_proteins(atoms, crystal_protein_atoms)
    if R is None:
        return None
    # Transform predicted ligand into crystal frame
    pred_conf = pred_lig_mol.GetConformer()
    new_coords = np.array([[pred_conf.GetAtomPosition(i).x,
                            pred_conf.GetAtomPosition(i).y,
                            pred_conf.GetAtomPosition(i).z]
                           for i in range(pred_lig_mol.GetNumAtoms())])
    aligned = transform_coords(new_coords, R, t)
    for i, (x, y, z) in enumerate(aligned):
        pred_conf.SetAtomPosition(i, (float(x), float(y), float(z)))
    rmsd = sym_corrected_rmsd(pred_lig_mol, crystal_lig_mol)
    if rmsd is None:
        return None
    return {"n_ca_aligned": n_ca, "lig_rmsd": rmsd, "lddt_pli_proxy": 1.0 / (1.0 + rmsd)}


def main(pred_root: Path):
    out_rows = []
    # Handle two layouts:
    #   1. Flat: <root>/<holo_id>/*.pdb
    #   2. Boltz: <root>/boltz_results_<holo_id>/predictions/<holo_id>/*.pdb
    candidates = sorted(d for d in pred_root.iterdir() if d.is_dir())
    holo_dirs = []
    for d in candidates:
        if d.name.startswith("boltz_results_"):
            holo_id = d.name[len("boltz_results_"):]
            inner = d / "predictions" / holo_id
            if inner.exists():
                holo_dirs.append((holo_id, inner))
        else:
            holo_dirs.append((d.name, d))
    print(f"Found {len(holo_dirs)} holo prediction dirs in {pred_root}")
    for holo_id, hd in holo_dirs:
        crystal_pdb = CRYSTAL / f"{holo_id}.pdb"
        if not crystal_pdb.exists():
            print(f"  {holo_id}: no crystal LIG file, skip")
            continue
        ctext = crystal_pdb.read_text()
        _, chetatms, cconects = parse_pdb_records(ctext)
        if not chetatms:
            print(f"  {holo_id}: no HETATM in crystal file")
            continue
        # Crystal protein atoms aren't in the LIG-only file — we need full holo
        # Look in data/external/pxr_holo, pxr_crystals, or pxr_crystals_chembl (nb16 expansion)
        holo_full = None
        for d in [REPO / "data" / "external" / "pxr_holo",
                  REPO / "data" / "external" / "pxr_crystals",
                  REPO / "data" / "external" / "pxr_crystals_chembl"]:
            cand = d / f"{holo_id}.pdb"
            if cand.exists():
                holo_full = cand
                break
        if holo_full is None:
            print(f"  {holo_id}: no full holo file found, skip")
            continue
        full_text = holo_full.read_text()
        crystal_atoms, _, _ = parse_pdb_records(full_text)
        # Build crystal lig mol from the LIG-only file's hetatms
        crystal_lig_mol = lig_from_records(chetatms, cconects)
        if crystal_lig_mol is None:
            print(f"  {holo_id}: crystal lig won't parse as mol")
            continue
        # Score every pred .pdb under this dir
        n_scored = 0
        for pdb in sorted(hd.rglob("*.pdb")):
            res = score_one(pdb, crystal_lig_mol, crystal_atoms)
            if res:
                try:
                    rel = str(pdb.relative_to(REPO))
                except ValueError:
                    rel = str(pdb)
                out_rows.append({
                    "holo_id": holo_id,
                    "pred_path": rel,
                    **res,
                })
                n_scored += 1
        print(f"  {holo_id}: {n_scored} preds scored, "
              f"best_rmsd={min((r['lig_rmsd'] for r in out_rows if r['holo_id'] == holo_id), default=None)}")
    if not out_rows:
        print("No scores written.")
        return
    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nWrote {len(out_rows)} scores to {OUT_CSV}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: score_validation.py PREDICTIONS_DIR (dir with one subdir per holo)")
        sys.exit(1)
    main(Path(sys.argv[1]).resolve())
