"""
Flexible N-model z-hybrid builder. Assemble a z-hybrid from ANY subset of pose sources:
  B2  OF3  AF3  PROT  CHAI  GNINA  GLIDE
Per model take each compound's best-pose confidence, z-score WITHIN model across the 184
(auto-calibrates scales), pick the max-z model's pose per compound. The 3-model (B2+OF3+AF3)
version scored 0.5414 (rank #2). Each independent/orthogonal source adds complementary poses.

Usage:
  python scripts/build_nmodel_zhybrid.py B2,OF3,AF3,PROT,GNINA  [out_tag]
  python scripts/build_nmodel_zhybrid.py B2,OF3,AF3,CHAI        chai5
Output: submissions/<out_tag>.zip + data/external/curation/<out_tag>_selection.csv
"""
from __future__ import annotations
import os, sys, json, glob, shutil, zipfile, csv
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COFOLD = os.path.join(REPO, "data", "external", "dargason_cofold")
B2_DIR = os.path.join(COFOLD, "predictions")
OF3_DIR = os.path.join(COFOLD, "openfold3_extract", "predictions")
N_PROTEIN = 293

# model -> best-pose dir for the generic (ranked_best) sources
RANKED_DIRS = {
    "AF3":   os.path.join(COFOLD, "af3_best"),
    "PROT":  os.path.join(COFOLD, "protenix_best"),
    "CHAI":  os.path.join(COFOLD, "chai_best"),
    # IntelliFold (IntFold) v2-flash best-of-5 (ranking = IntFold ranking_score). GT-tested
    # to beat the cofolder pool on 5/8 PXR failure holos; folded in as a failure-tail rescuer.
    "INTFOLD": os.path.join(COFOLD, "intfold_best"),
    "GNINA": os.path.join(COFOLD, "gnina_best"),
    "GLIDE": os.path.join(COFOLD, "glide_best"),
    # Boltz-2.1 hosted-API best-of-5 (ranking = -complex_ipde). restr = per-ligand dossier pocket
    # restraints; free = unrestrained. Both produced by scripts/boltz_run_184.py extract.
    "B2API":      os.path.join(REPO, "data", "external", "boltz_api", "run184_restr_best"),
    "B2API_FREE": os.path.join(REPO, "data", "external", "boltz_api", "run184_free_best"),
}


def b2_best(sid):
    best = None
    for pdb in glob.glob(os.path.join(B2_DIR, sid, "**", "*.pdb"), recursive=True):
        c = os.path.join(os.path.dirname(pdb), f"confidence_{os.path.basename(pdb)[:-4]}.json")
        if not os.path.exists(c): continue
        try: v = json.load(open(c)).get("complex_ipde")
        except: continue
        if v is None: continue
        s = -float(v)
        if best is None or s > best[1]: best = (pdb, s)
    return best


def of3_best(sid):
    best = None
    for pdb in glob.glob(os.path.join(OF3_DIR, sid, "**", "*.pdb"), recursive=True):
        stem = os.path.basename(pdb)[:-4]; base = stem.replace("_model", ""); d = os.path.dirname(pdb)
        full = os.path.join(d, f"{base}_confidences.json"); agg = os.path.join(d, f"{base}_confidences_aggregated.json")
        s = None
        try:
            if os.path.exists(full):
                pae = np.array(json.load(open(full))["pae"]); n = pae.shape[0]
                if n - N_PROTEIN > 0: s = float(-pae[:N_PROTEIN, N_PROTEIN:n].mean())
            elif os.path.exists(agg):
                iptm = json.load(open(agg)).get("chain_pair_iptm", {}).get("(A, L)")
                if iptm is not None: s = float(iptm)
        except: s = None
        if s is None: continue
        if best is None or s > best[1]: best = (pdb, s)
    return best


_rank_cache = {}
def ranked_best(model, sid):
    d = RANKED_DIRS[model]
    if d not in _rank_cache:
        rj = os.path.join(d, "_ranking.json")
        _rank_cache[d] = json.load(open(rj)) if os.path.exists(rj) else {}
    rk = _rank_cache[d]
    p = os.path.join(d, f"{sid}.pdb")
    if os.path.exists(p) and rk.get(sid) is not None:
        return (p, float(rk[sid]))
    return None


def get_best(model, sid):
    if model == "B2": return b2_best(sid)
    if model == "OF3": return of3_best(sid)
    return ranked_best(model, sid)


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    models = [m.strip().upper() for m in sys.argv[1].split(",") if m.strip()]
    tag = sys.argv[2] if len(sys.argv) > 2 else "zhybrid_" + "_".join(m.lower() for m in models)
    out_dir = os.path.join(REPO, "submissions", tag)
    out_zip = os.path.join(REPO, "submissions", tag + ".zip")
    sids = sorted(d for d in os.listdir(B2_DIR) if os.path.isdir(os.path.join(B2_DIR, d)))

    M = {m: {} for m in models}
    for sid in sids:
        for m in models:
            r = get_best(m, sid)
            if r: M[m][sid] = r
    print("poses available per model:", {m: len(M[m]) for m in models})

    def zmap(d):
        if not d: return {}
        v = np.array([d[s][1] for s in d]); mu, sd = v.mean(), v.std() or 1.0
        return {s: (d[s][1] - mu) / sd for s in d}
    Z = {m: zmap(M[m]) for m in models}

    if os.path.exists(out_dir): shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    wins = {m: 0 for m in models}; rows = []
    for sid in sids:
        cands = [(m, Z[m][sid], M[m][sid][0]) for m in models if sid in Z[m]]
        if not cands: continue
        m, z, path = max(cands, key=lambda c: c[1])
        shutil.copy(path, os.path.join(out_dir, f"{sid}.pdb"))
        wins[m] += 1; rows.append((sid, m, round(float(z), 2)))
    print("model wins:", wins, " total:", sum(wins.values()))
    with open(os.path.join(REPO, "data", "external", "curation", tag + "_selection.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["structure", "model", "z"]); w.writerows(rows)
    pdbs = sorted(glob.glob(os.path.join(out_dir, "*.pdb")))
    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in pdbs: zf.write(p, arcname=os.path.basename(p))
    print(f"zipped {len(pdbs)} -> {out_zip} ({os.path.getsize(out_zip)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
