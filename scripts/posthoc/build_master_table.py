#!/usr/bin/env python3
"""
Post-hoc: master cross-model per-pose table for all 184 PXR test ligands.

For each ligand, load every available model's best pose, put them all in a common
protein frame (CA superposition), and quantify:
  - cross-model ligand-pose disagreement (mean pairwise RMSD)  -> difficulty proxy
  - consensus medoid model + per-model deviation from medoid
  - per-model confidence surrogate (mean ligand B-factor = pLDDT)
Merged with ligand properties (MW, rotB, fragment/analog, series, pEC50).

No ground truth exists for the 184; cross-model agreement is the difficulty signal.
Output: data/processed/posthoc/master_184.parquet + per-model long table.
"""
import os
import sys
import glob
import itertools
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd

sys.path.insert(0, "scripts")
import pose_lib as pl

BASE = "data/external"
MODELS = {
    "boltz1": "boltz1_184_best",
    "of3ft": "of3ft_best",
    "protenix_v2": "openprotein_protenix_v2_best",
    "decaf": "decaf_best",
    "esm2_full": "esm2_full_best",
    "esm2_shallow": "esm2_shallow_best",
    "esm2_none": "esm2_none_best",
    "boltz1x": "openprotein_boltz1x_best",
    "esmfold2fast": "openprotein_esmfold2fast_best",
    "rf3": "openprotein_rf3_best",
    "apo": "apo_best",
    # partial-coverage models (included where present)
    "af3_iface": "af3_iface_best",
    "af3_massive": "af3_massive_best",
    "chai": "chai_massive_best",
    "restr": "restr_best",
}
REF = "boltz1"  # reference frame


def elem_matched_rmsd(xyzA, elemsA, xyzB, elemsB):
    """Greedy nearest-neighbour RMSD within matching elements (symmetry-agnostic).
    Adequate for cross-model disagreement; not a substitute for GT RMSD."""
    if len(elemsA) != len(elemsB):
        # match by element multiset on the smaller set
        pass
    used = set()
    sq = []
    for i, (a, ea) in enumerate(zip(xyzA, elemsA)):
        best = None
        bj = -1
        for j, (b, eb) in enumerate(zip(xyzB, elemsB)):
            if j in used or eb != ea:
                continue
            d = float(np.sum((a - b) ** 2))
            if best is None or d < best:
                best = d
                bj = j
        if bj >= 0:
            used.add(bj)
            sq.append(best)
    if not sq:
        return np.nan
    return float(np.sqrt(np.mean(sq)))


_SKIP_HET = {"HOH", "WAT", "NA", "CL", "K", "MG", "ZN", "CA", "SO4", "PO4", "GOL", "EDO", "ACT"}


def parse_pose_robust(path):
    """Robust to ligand convention (LIG/chain B, l01/chain X, L:0/chain X, ...).
    CA from ATOM CA lines; ligand = any non-water/ion HETATM."""
    ca, ca_res = {}, {}
    lig_xyz, lig_elems, lig_bfac = [], [], []
    for ln in open(path):
        rec = ln[:6].strip()
        if rec == "ATOM" and ln[12:16].strip() == "CA":
            try:
                rs = int(ln[22:26]); xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                ca[rs] = xyz; ca_res[rs] = ln[17:20].strip()
            except Exception:
                pass
        elif rec == "HETATM":
            resn = ln[17:20].strip().upper()
            if resn in _SKIP_HET:
                continue
            try:
                xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                el = ln[76:78].strip() or "".join(c for c in ln[12:16].strip() if c.isalpha())[:1]
                bf = float(ln[60:66]) if ln[60:66].strip() else np.nan
                lig_xyz.append(xyz); lig_elems.append(el.capitalize()); lig_bfac.append(bf)
            except Exception:
                pass
    import numpy as _np
    return {"ca": {k: _np.array(v) for k, v in ca.items()}, "ca_res": ca_res,
            "lig_xyz": _np.array(lig_xyz), "lig_elems": lig_elems,
            "lig_bfac": _np.array(lig_bfac) if lig_bfac else _np.array([_np.nan])}


def load_poses(sid):
    poses = {}
    for m, d in MODELS.items():
        f = os.path.join(BASE, d, f"{sid}.pdb")
        if os.path.exists(f) and os.path.getsize(f) > 1000:
            try:
                p = parse_pose_robust(f)
                if len(p["lig_xyz"]) >= 4 and len(p["ca"]) >= 50:
                    poses[m] = p
            except Exception:
                pass
    return poses


def align_to_ref(pose, ref_pose):
    """Return ref-frame ligand xyz for `pose`, aligning its CA onto ref CA.
    align_ca returns (R, t, rms, n_keys); align_by_sequence used as numbering fallback."""
    res = pl.align_ca(pose["ca"], ref_pose["ca"])
    if res is not None:
        R, t = res[0], res[1]
        return pl.apply_T(R, t, pose["lig_xyz"])
    try:
        r2 = pl.align_by_sequence(pose["ca"], pose["ca_res"], ref_pose["ca"], ref_pose["ca_res"])
        R, t = r2[0], r2[1]
        return pl.apply_T(R, t, pose["lig_xyz"])
    except Exception:
        return None


def main():
    props = pd.read_csv("data/external/curation/master_ligands.csv")
    sids = props["structure"].tolist()
    os.makedirs("data/processed/posthoc", exist_ok=True)

    rows = []       # per-ligand summary
    long_rows = []  # per-ligand-per-model
    for n, sid in enumerate(sids):
        poses = load_poses(sid)
        if REF not in poses or len(poses) < 3:
            rows.append({"structure": sid, "n_models": len(poses)})
            continue
        ref = poses[REF]
        # ref-frame ligand coords per model
        frame = {}
        for m, p in poses.items():
            xyz = ref["lig_xyz"] if m == REF else align_to_ref(p, ref)
            if xyz is not None and len(xyz) == len(p["lig_elems"]):
                frame[m] = (xyz, p["lig_elems"], float(np.nanmean(p["lig_bfac"])))
        models = list(frame)
        # pairwise RMSD matrix
        pw = {}
        for a, b in itertools.combinations(models, 2):
            r = elem_matched_rmsd(frame[a][0], frame[a][1], frame[b][0], frame[b][1])
            pw[(a, b)] = pw[(b, a)] = r
        # medoid: min mean RMSD to others
        mean_to_others = {}
        for a in models:
            vals = [pw[(a, b)] for b in models if b != a and not np.isnan(pw[(a, b)])]
            mean_to_others[a] = np.mean(vals) if vals else np.nan
        medoid = min(mean_to_others, key=lambda k: (np.isnan(mean_to_others[k]), mean_to_others[k]))
        allpw = [v for k, v in pw.items() if k[0] < k[1] and not np.isnan(v)]
        disagreement = np.mean(allpw) if allpw else np.nan
        # count models within 2A of medoid (consensus cluster size)
        consensus_n = sum(1 for a in models if a != medoid and not np.isnan(pw[(a, medoid)]) and pw[(a, medoid)] <= 2.0) + 1
        rows.append({
            "structure": sid,
            "n_models": len(frame),
            "disagreement_mean_pw_rmsd": round(disagreement, 3) if not np.isnan(disagreement) else np.nan,
            "disagreement_max_pw_rmsd": round(np.nanmax(allpw), 3) if allpw else np.nan,
            "medoid_model": medoid,
            "consensus_cluster_n": consensus_n,
            "consensus_frac": round(consensus_n / len(frame), 3),
            "mean_lig_plddt": round(np.mean([frame[m][2] for m in models]), 1),
        })
        for m in models:
            if m == medoid:
                rtm = 0.0
            elif not np.isnan(pw[(m, medoid)]):
                rtm = round(pw[(m, medoid)], 3)
            else:
                rtm = np.nan  # failed comparison -> NaN, not 0 (bugfix)
            long_rows.append({
                "structure": sid, "model": m, "rmsd_to_medoid": rtm,
                "lig_plddt": round(frame[m][2], 1) if not np.isnan(frame[m][2]) else np.nan,
                "is_medoid": m == medoid,
            })
        if n % 30 == 0:
            print(f"  {n}/184 done (last {sid}: {len(frame)} models, disagree={disagreement:.2f})")

    summ = pd.DataFrame(rows).merge(props, on="structure", how="left")
    summ.to_parquet("data/processed/posthoc/master_184.parquet")
    summ.to_csv("data/processed/posthoc/master_184.csv", index=False)
    pd.DataFrame(long_rows).to_parquet("data/processed/posthoc/master_184_long.parquet")
    ok = summ["disagreement_mean_pw_rmsd"].notna()
    print(f"\nDONE. {ok.sum()}/184 ligands with cross-model geometry.")
    print(f"disagreement (mean pairwise RMSD): median={summ.loc[ok,'disagreement_mean_pw_rmsd'].median():.2f}A")
    print(f"medoid model distribution:\n{summ.loc[ok,'medoid_model'].value_counts().to_string()}")
    print(f"consensus_frac by population:")
    print(summ[ok].groupby('population')['consensus_frac'].median().to_string())


if __name__ == "__main__":
    main()
