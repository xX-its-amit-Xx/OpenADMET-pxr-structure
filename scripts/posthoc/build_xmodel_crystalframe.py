#!/usr/bin/env python3
"""
Correct cross-model geometry in a CRYSTAL reference frame (all 184).

Why: the original master table aligned every model to `boltz1`, but boltz1's 184
protein exports fold ~22 A from the real PXR LBD (systematic: 0/12 within 3 A),
along with decaf. Aligning to a misfolded reference inflated the disagreement
metric ~2x and made "Boltz is an outlier" partly an alignment artefact.

Fix: use a real PXR LBD crystal (1ILH) as the universal frame. Align each model's
CA onto it (residue-identity offset); keep only well-folded models (< 3 A). Recompute
disagreement / medoid / consensus / rmsd_to_medoid from element-matched ligand RMSD in
that frame. Also emit per-compound aligned ligand poses for the 3D scatter, plus the
shared crystal pocket backdrop.

Out:
  data/processed/posthoc/xmodel_crystalframe.parquet   (per-compound corrected summary)
  data/processed/posthoc/xmodel_long.parquet           (per-model rmsd_to_medoid, folded_ok)
  data/processed/posthoc/journey3d/<sid>.json          (aligned ligand poses, 3D)
  docs/journey_pocket.pdb                               (shared crystal pocket + ligand backdrop)
"""
import os, sys, json, glob, itertools
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, "scripts")
import pose_lib as pl

HEAVY = set("C N O P S F Cl Br I B Se".split())
MODELS = {'boltz1':'boltz1_184_best','of3ft':'of3ft_best','protenix_v2':'openprotein_protenix_v2_best',
 'decaf':'decaf_best','esm2_full':'esm2_full_best','esm2_shallow':'esm2_shallow_best','esm2_none':'esm2_none_best',
 'boltz1x':'openprotein_boltz1x_best','esmfold2fast':'openprotein_esmfold2fast_best','rf3':'openprotein_rf3_best',
 'apo':'apo_best','af3_iface':'af3_iface_best','af3_massive':'af3_massive_best','chai':'chai_massive_best','restr':'restr_best'}
BASE = "data/external"
CRYSTAL = "data/external/pxr_holo/1ilh.pdb"; LIGCODE = "SRL"
FOLD_TOL = 3.0
OUT3D = "data/processed/posthoc/journey3d"
ML = pd.read_parquet("data/processed/posthoc/master_184_long.parquet")


def emrmsd(a, ae, b, be):
    if len(ae) != len(be) or sorted(ae) != sorted(be):
        return np.nan
    D = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    big = D.max() + 1e4
    C = D.copy(); C[np.array(ae)[:, None] != np.array(be)[None, :]] = big
    r, c = linear_sum_assignment(C)
    return float(np.sqrt((D[r, c] ** 2).mean()))


def main():
    os.makedirs(OUT3D, exist_ok=True)
    cr = pl.parse_crystal(CRYSTAL, LIGCODE)
    ch = max(cr["chains"], key=lambda c: len(cr["chains"][c]))
    cca, cres = cr["chains"][ch], cr["chains_res"][ch]
    clig = np.array(cr["lig_copies"][0]["xyz"]); clc = clig.mean(0)
    # shared backdrop: pocket residues (CA+sidechain of residues within 12A of crystal ligand) + crystal ligand
    pocket_keys = {k for k in cca if np.linalg.norm(np.array(cca[k]) - clc) < 12}
    pocket_lines = [ln.rstrip("\n") for ln in open(CRYSTAL)
                    if ln[:4] == "ATOM" and ln[21] == ch and int(ln[22:26]) in pocket_keys]
    crystal_lig_lines = [ln.rstrip("\n") for ln in open(CRYSTAL)
                         if ln[:6] == "HETATM" and ln[17:20].strip() == LIGCODE]
    open("docs/journey_pocket.pdb", "w").write("\n".join(pocket_lines) + "\n")
    open("docs/journey_crystal_lig.pdb", "w").write("\n".join(crystal_lig_lines) + "\n")
    print(f"backdrop: {len(pocket_keys)} pocket residues, {len(crystal_lig_lines)} crystal ligand atoms", flush=True)

    lig_plddt = {(r.structure, r.model): float(r.lig_plddt) for r in ML.itertuples()}
    sids = pd.read_parquet("data/processed/posthoc/master_184.parquet")["structure"].tolist()

    summ, longrows = [], []
    for n, sid in enumerate(sids):
        frame = {}          # model -> (ligT[heavy], elems[heavy])
        excluded = []
        for m, d in MODELS.items():
            p = os.path.join(BASE, d, f"{sid}.pdb")
            if not os.path.exists(p):
                continue
            try:
                pp = pl.parse_pose(p)
            except Exception:
                continue
            if len(pp["ca"]) < 60 or len(pp["lig_elems"]) == 0:
                continue
            al = pl.align_by_sequence(pp["ca"], pp["ca_res"], cca, cres)
            if al is None or al[2] > FOLD_TOL:
                excluded.append(m); continue
            R, t, prms, _ = al
            ligT = pl.apply_T(R, t, pp["lig_xyz"])
            keep = [i for i, e in enumerate(pp["lig_elems"]) if e in HEAVY]
            frame[m] = (ligT[keep], [pp["lig_elems"][i] for i in keep], prms)
        models = list(frame)
        if len(models) < 3:
            summ.append(dict(structure=sid, n_wellfolded=len(models), excluded=";".join(excluded),
                             disagreement_cf=np.nan, medoid_cf=None, consensus_frac_cf=np.nan)); continue
        # pairwise element-matched RMSD
        pw = {}
        for a, b in itertools.combinations(models, 2):
            r = emrmsd(frame[a][0], frame[a][1], frame[b][0], frame[b][1])
            pw[(a, b)] = pw[(b, a)] = r
        mean_to = {a: np.nanmean([pw[(a, b)] for b in models if b != a]) for a in models}
        medoid = min(mean_to, key=lambda k: (np.isnan(mean_to[k]), mean_to[k]))
        allpw = [v for k, v in pw.items() if k[0] < k[1] and not np.isnan(v)]
        disagreement = float(np.mean(allpw)) if allpw else np.nan
        cons = (sum(1 for a in models if a != medoid and not np.isnan(pw[(a, medoid)]) and pw[(a, medoid)] <= 2.0) + 1) / len(models)
        # plddt pick among well-folded
        pk = max(models, key=lambda m: lig_plddt.get((sid, m), -1))
        for m in models:
            rmd = float(np.nanmean([pw[(m, b)] for b in models if b != m])) if len(models) > 1 else 0.0
            longrows.append(dict(structure=sid, model=m, rmsd_to_medoid_cf=round(float(pw[(m, medoid)] if m != medoid else 0.0), 2),
                                 mean_pw_cf=round(rmd, 2), lig_plddt=lig_plddt.get((sid, m)),
                                 is_medoid=(m == medoid), is_pick=(m == pk), prot_rmsd=round(frame[m][2], 2)))
        summ.append(dict(structure=sid, n_wellfolded=len(models), excluded=";".join(excluded),
                         disagreement_cf=round(disagreement, 2), medoid_cf=medoid,
                         consensus_frac_cf=round(cons, 2), plddt_pick_cf=pk,
                         selection_decoupled_cf=bool(pk != medoid)))
        # ---- 3D scatter json (well-folded models, crystal frame) ----
        j = dict(sid=sid, medoid=medoid, plddt_pick=pk, excluded=excluded,
                 models=[dict(model=m, lig=[[round(float(v), 2) for v in a] for a in frame[m][0].tolist()],
                              el=frame[m][1], rmsd_to_medoid=round(float(pw[(m, medoid)] if m != medoid else 0.0), 2),
                              lig_plddt=lig_plddt.get((sid, m)), is_medoid=(m == medoid), is_pick=(m == pk),
                              prot_rmsd=round(frame[m][2], 2)) for m in models])
        j["models"].sort(key=lambda x: x["rmsd_to_medoid"])
        json.dump(j, open(os.path.join(OUT3D, f"{sid}.json"), "w"))
        if n % 20 == 0:
            print(f"  {n+1}/{len(sids)} {sid}: {len(models)} folded, {len(excluded)} excluded, disag {disagreement:.2f}", flush=True)
    pd.DataFrame(summ).to_parquet("data/processed/posthoc/xmodel_crystalframe.parquet")
    pd.DataFrame(longrows).to_parquet("data/processed/posthoc/xmodel_long.parquet")
    # provenance
    S = pd.DataFrame(summ)
    print(f"\nDONE. {len(S)} compounds | median n_wellfolded {S.n_wellfolded.median():.0f}")
    print(f"median disagreement crystal-frame {S.disagreement_cf.median():.2f} A (was 4.86 in boltz1-frame)")
    from collections import Counter
    exc = Counter(x for row in summ for x in (row["excluded"].split(";") if row["excluded"] else []))
    print("most-excluded models:", dict(exc.most_common(6)))


if __name__ == "__main__":
    main()
