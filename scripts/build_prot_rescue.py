"""Parameterized Protenix failure-tail rescue builder.

Base = 4-model z-hybrid (AF3,B2,OF3,CHAI) at submissions/ai_4model_base/.
Rank the 184 by POOL confidence = max within-model z across the 4 base models (the same
score the base z-hybrid selects on). The N LOWEST (the failure tail where the pool is least
sure of its pose) get their pose SWAPPED to the Protenix-v2 best pose; all others keep base.

This reproduces the proven targeted-rescue family:
  N=12 -> prot_rescue   = 0.5629 (NEW BEST rank 2)
  N=20 -> prot_rescue20 = 0.5587 (regressed -> optimum is <=12)
Use this to bracket the optimal swap count (N=8, N=10, ...).

Usage: python scripts/build_prot_rescue.py N TAG
  e.g.  python scripts/build_prot_rescue.py 8 prot_rescue8
"""
from __future__ import annotations
import os, sys, glob, shutil, zipfile, csv
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import build_nmodel_zhybrid as B

BASE_DIR = os.path.join(REPO, "submissions", "ai_4model_base")
# PROT_DIR override (env) lets us swap in the best-of-100 high-sampling failtail poses (data/external/failtail_best)
PROT_DIR = os.environ.get("PROT_DIR", os.path.join(REPO, "data", "external", "dargason_cofold", "protenix_best"))
MODELS = ["AF3", "B2", "OF3", "CHAI"]


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    N = int(sys.argv[1]); tag = sys.argv[2]
    assert os.path.isdir(BASE_DIR), f"missing base {BASE_DIR}"

    sids = sorted(d for d in os.listdir(B.B2_DIR) if os.path.isdir(os.path.join(B.B2_DIR, d)))
    # per-model within-model z of best-pose confidence (same calibration the base uses)
    M = {m: {} for m in MODELS}
    for sid in sids:
        for m in MODELS:
            r = B.get_best(m, sid)
            if r: M[m][sid] = r

    def zmap(d):
        if not d: return {}
        v = np.array([d[s][1] for s in d]); mu, sd = v.mean(), v.std() or 1.0
        return {s: (d[s][1] - mu) / sd for s in d}
    Z = {m: zmap(M[m]) for m in MODELS}

    # pool confidence per ligand = max within-model z across the 4 base models
    pool = {}
    for sid in sids:
        zs = [Z[m][sid] for m in MODELS if sid in Z[m]]
        if zs: pool[sid] = max(zs)

    # N lowest pool confidence = the failure tail to rescue with Protenix
    ranked = sorted(pool, key=lambda s: pool[s])
    swap = [s for s in ranked if os.path.isfile(os.path.join(PROT_DIR, f"{s}.pdb"))
            and os.path.isfile(os.path.join(BASE_DIR, f"{s}.pdb"))][:N]
    swap_set = set(swap)

    out_dir = os.path.join(REPO, "submissions", tag)
    if os.path.exists(out_dir): shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    nb = 0
    for sid in sids:
        src = os.path.join(PROT_DIR if sid in swap_set else BASE_DIR, f"{sid}.pdb")
        if not os.path.isfile(src):
            src = os.path.join(BASE_DIR, f"{sid}.pdb")
        if os.path.isfile(src):
            shutil.copy(src, os.path.join(out_dir, f"{sid}.pdb")); nb += 1
    print(f"N={N} swapped {len(swap)} to Protenix (lowest pool-z): {swap}")
    with open(os.path.join(REPO, "data", "external", "curation", tag + "_selection.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["structure", "source", "pool_z"])
        for sid in sids:
            w.writerow([sid, "PROT" if sid in swap_set else "BASE", round(pool.get(sid, 0.0), 3)])
    pdbs = sorted(glob.glob(os.path.join(out_dir, "*.pdb")))
    out_zip = os.path.join(REPO, "submissions", tag + ".zip")
    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in pdbs: zf.write(p, arcname=os.path.basename(p))
    print(f"zipped {len(pdbs)} -> {out_zip} ({os.path.getsize(out_zip)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
