"""
3-model z-hybrid selection: Boltz-2 + OpenFold3 + AF3.
Per model, take each compound's best-pose confidence, z-score WITHIN model across the 184
(this auto-calibrates the different confidence scales -- the mechanism behind #262/#13's
B2+OF3 z-hybrid that scored 0.4997). Per compound, output the pose of the model with the
highest z. AF3 best pose + ranking_score come from the Modal run (af3_best/).

Output: submissions/ai_4model_zhybrid/  + .zip   (HOLD -- validate before submitting)
"""
from __future__ import annotations
import os, sys, json, glob, shutil, zipfile, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COFOLD = os.path.join(REPO, "data", "external", "dargason_cofold")
B2_DIR = os.path.join(COFOLD, "predictions")
OF3_DIR = os.path.join(COFOLD, "openfold3_extract", "predictions")
AF3_BEST = os.path.join(COFOLD, "af3_best")
PROT_BEST = os.path.join(COFOLD, "protenix_best")
N_PROTEIN = 293
OUT_DIR = os.path.join(REPO, "submissions", "ai_4model_zhybrid")
OUT_ZIP = os.path.join(REPO, "submissions", "ai_4model_zhybrid.zip")

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

def main():
    sids = sorted(d for d in os.listdir(B2_DIR) if os.path.isdir(os.path.join(B2_DIR, d)))
    af3_rank = json.load(open(os.path.join(AF3_BEST, "_ranking.json")))
    # gather per-model best (path, signal)
    B, O, A, P = {}, {}, {}, {}
    prot_rank = json.load(open(os.path.join(PROT_BEST, "_ranking.json"))) if os.path.exists(os.path.join(PROT_BEST, "_ranking.json")) else {}
    for sid in sids:
        b = b2_best(sid); o = of3_best(sid)
        if b: B[sid] = b
        if o: O[sid] = o
        ap = os.path.join(AF3_BEST, f"{sid}.pdb")
        if os.path.exists(ap) and af3_rank.get(sid) is not None:
            A[sid] = (ap, float(af3_rank[sid]))
        pp = os.path.join(PROT_BEST, f"{sid}.pdb")
        if os.path.exists(pp) and prot_rank.get(sid) is not None:
            P[sid] = (pp, float(prot_rank[sid]))
    # within-model z
    def zmap(M):
        v = np.array([M[s][1] for s in M]); mu, sd = v.mean(), v.std() or 1.0
        return {s: (M[s][1] - mu) / sd for s in M}
    zB, zO, zA, zP = zmap(B), zmap(O), zmap(A), zmap(P)
    if os.path.exists(OUT_DIR): shutil.rmtree(OUT_DIR)
    os.makedirs(OUT_DIR)
    wins = {"B2": 0, "OF3": 0, "AF3": 0, "PROT": 0}; rows = []
    for sid in sids:
        cands = []
        if sid in zB: cands.append(("B2", zB[sid], B[sid][0]))
        if sid in zO: cands.append(("OF3", zO[sid], O[sid][0]))
        if sid in zA: cands.append(("AF3", zA[sid], A[sid][0]))
        if sid in zP: cands.append(("PROT", zP[sid], P[sid][0]))
        if not cands: continue
        m, z, path = max(cands, key=lambda c: c[1])
        shutil.copy(path, os.path.join(OUT_DIR, f"{sid}.pdb"))
        wins[m] += 1; rows.append((sid, m, round(z, 2)))
    print("model wins:", wins, " total:", sum(wins.values()))
    with open(os.path.join(REPO, "data", "external", "curation", "zhybrid4_selection.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["structure", "model", "z"]); w.writerows(rows)
    pdbs = sorted(glob.glob(os.path.join(OUT_DIR, "*.pdb")))
    with zipfile.ZipFile(OUT_ZIP, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in pdbs: zf.write(p, arcname=os.path.basename(p))
    print(f"zipped {len(pdbs)} -> {OUT_ZIP} ({os.path.getsize(OUT_ZIP)/1e6:.2f} MB)")
    # divergence vs #262 (all-B2/OF3 IPDE); here just report AF3 incursion
    print(f"AF3 won {wins['AF3']}/{sum(wins.values())} compounds (these differ from B2+OF3-only #262)")

if __name__ == "__main__":
    main()
