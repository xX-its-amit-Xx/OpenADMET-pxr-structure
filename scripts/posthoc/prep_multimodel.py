#!/usr/bin/env python3
"""
Multi-model divergence scene: for one holo, put every architecture's co-folded
ligand into the true crystal frame and emit JSON for the 3Dmol viewer. Shows the
user's requested "different models -> different poses" against real crystal GT,
and makes the input-fidelity finding visible (all architectures miss under a
generic-sequence setup). Reuses the authoritative scoring path from score_cofold.
"""
import os, sys, glob, json
import numpy as np
import pandas as pd

sys.path.insert(0, "scripts"); sys.path.insert(0, "scripts/posthoc")
import pose_lib as pl
from score_cofold import (expected_lig, heavy_only, split_copies, elem_rmsd,
                          crystal_lig_atoms, crystal_path)

HOLO = os.environ.get("MM_HOLO", "1M13")
COFOLD = "data/external/posthoc_holo_cofold"


def main():
    hr = pd.read_csv("data/processed/validation_set/holos_fixed.csv")
    r = hr[hr.holo_id == HOLO].iloc[0]
    lig_code = r["lig_code"]
    exp = expected_lig(r.get("smiles_canonical") or r.get("smiles"))
    n_heavy, ref_el = exp
    cpath = crystal_path(HOLO)
    cr = pl.parse_crystal(cpath, lig_code)
    chain = max(cr["chains"], key=lambda c: len(cr["chains"][c]))
    cca, cres = cr["chains"][chain], cr["chains_res"][chain]
    cxyz, cel = crystal_lig_atoms(cpath, lig_code)

    # crystal protein ATOM block + a single crystal ligand copy (first n_heavy heavy atoms
    # of the copy nearest the pocket centroid) for display
    prot_lines = [ln.rstrip("\n") for ln in open(cpath) if ln[:4] == "ATOM" and ln[21] == chain]
    # crystal ligand: take one representative copy (first n_heavy atoms matching ref)
    cryst_lig = {"xyz": cxyz[:n_heavy].round(3).tolist(), "el": cel[:n_heavy]}

    models = []
    for model in sorted(os.listdir(COFOLD)):
        pf = glob.glob(f"{COFOLD}/{model}/{HOLO}_s0.pdb")
        if not pf:
            continue
        p = pl.parse_pose(pf[0])
        al = pl.align_by_sequence(p["ca"], p["ca_res"], cca, cres)
        if al is None or al[2] > 3.0:
            continue
        R, t, prms, _ = al
        hx, he = heavy_only(p["lig_xyz"], list(p["lig_elems"]))
        samples = split_copies(hx, he, n_heavy, ref_el)
        if not samples:
            continue
        sx, se = samples[0]
        sxT = pl.apply_T(R, t, sx)
        rmsd = elem_rmsd(sxT, se, cxyz, cel)
        models.append({"model": model, "rmsd": round(rmsd, 2), "prot_rmsd": round(prms, 2),
                       "lig": sxT.round(3).tolist(), "el": se})
    models.sort(key=lambda m: m["rmsd"])
    out = {"holo": HOLO, "ligcode": lig_code, "n_heavy": n_heavy,
           "protein_pdb": "\n".join(prot_lines), "crystal_lig": cryst_lig,
           "models": models}
    os.makedirs("data/processed/posthoc", exist_ok=True)
    p = f"data/processed/posthoc/multimodel_{HOLO}.json"
    json.dump(out, open(p, "w"))
    print(f"wrote {p}: {len(models)} models on {HOLO} ({lig_code}, {n_heavy} heavy)")
    for m in models:
        print(f"  {m['model']:12s} lig-RMSD {m['rmsd']:5.2f}A  (protein {m['prot_rmsd']:.2f}A)")


if __name__ == "__main__":
    main()
