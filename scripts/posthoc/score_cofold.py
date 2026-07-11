#!/usr/bin/env python3
"""
Score the cross-model holo co-fold vs crystal ligands -> per-model TRUE RMSD.

Robustness: OpenProtein merged the N diffusion samples into ONE ligand block
(chain X, resSeq 1); crystals may contain multiple ligand copies. We compute the
expected heavy-atom count per ligand from the canonical SMILES (largest fragment),
split BOTH sides into copies of that size (heavy atoms only, hydrogens dropped),
align the co-fold protein onto the crystal by residue-identity sequence offset,
then element-matched RMSD (Hungarian) for every (co-fold sample x crystal copy),
taking the min over crystal copies. Oracle = best sample; also a plddt-style pick.

Outputs data/processed/posthoc/cofold_scores.csv (one row per model/holo/sample)
and a per-model summary printed to stdout.
"""
import os, sys, glob, json
import numpy as np
import pandas as pd
from collections import Counter
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, "scripts")
import pose_lib as pl

COFOLD = "data/external/posthoc_holo_cofold"
HOLOS = pd.read_csv("data/processed/validation_set/holos_fixed.csv")
CRYSTAL_ROOTS = ["data/external/pxr_holo", "data/external/rcsb_holo",
                 "data/external/novel_val/crystals", "data/external"]
HEAVY = set("C N O P S F Cl Br I B Se".split())
_CRYSTAL_CACHE = {}


def expected_lig(smiles):
    """(n_heavy, element Counter) for the largest fragment of a (possibly multi-component) SMILES."""
    from rdkit import Chem
    if not isinstance(smiles, str) or not smiles:
        return None
    frags = smiles.split(".")
    best = None
    for fr in frags:
        m = Chem.MolFromSmiles(fr)
        if m is None:
            continue
        el = Counter(a.GetSymbol() for a in m.GetAtoms())  # heavy only (implicit H)
        n = sum(el.values())
        if best is None or n > best[0]:
            best = (n, el)
    return best


def heavy_only(xyz, elems):
    keep = [i for i, e in enumerate(elems) if e in HEAVY]
    return xyz[keep], [elems[i] for i in keep]


def _elem_from_pdb_line(ln):
    e = ln[76:78].strip()
    if not e:
        e = ln[12:16].strip()[:1]
    return e[0].upper() + e[1:].lower() if len(e) > 1 else e.upper()


def crystal_lig_atoms(path, lig_code):
    """All heavy ligand atoms from the crystal (across all copies/residues). Element-
    constrained rectangular assignment against a single pose copy naturally scores to
    the nearest crystal copy, so no per-copy splitting is needed."""
    xyz, el = [], []
    for ln in open(path):
        if ln[:6] not in ("HETATM", "ATOM  "):
            continue
        if ln[17:20].strip() != lig_code:
            continue
        e = _elem_from_pdb_line(ln)
        if e not in HEAVY:
            continue
        xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        el.append(e)
    return np.asarray(xyz), el


def split_copies(xyz, elems, n_per, ref_el):
    """Split a merged atom block into contiguous copies of size n_per whose element
    multiset matches ref_el. Returns list of (xyz_k, elems_k) or [] if not clean."""
    if n_per == 0 or len(elems) % n_per != 0:
        return []
    k = len(elems) // n_per
    out = []
    for i in range(k):
        s = slice(i * n_per, (i + 1) * n_per)
        ek = elems[s]
        if Counter(ek) != ref_el:
            return []          # atom order differs across copies -> bail, caller logs
        out.append((xyz[s], ek))
    return out


def elem_rmsd(a_xyz, a_el, b_xyz, b_el):
    """Element-constrained min-RMSD of pose atoms a (small) onto crystal atoms b
    (>= a; may hold several copies). Rectangular Hungarian -> nearest-copy match."""
    D = np.linalg.norm(a_xyz[:, None, :] - b_xyz[None, :, :], axis=2)
    big = D.max() + 1e4
    C = D.copy()
    C[np.array(a_el)[:, None] != np.array(b_el)[None, :]] = big
    r, c = linear_sum_assignment(C)
    return float(np.sqrt((D[r, c] ** 2).mean()))


def crystal_path(h):
    if h in _CRYSTAL_CACHE:
        return _CRYSTAL_CACHE[h]
    hits = []
    for r in CRYSTAL_ROOTS:
        for name in (h.lower(), h.upper()):
            hits += glob.glob(f"{r}/**/{name}.pdb", recursive=True)
    hits = [x for x in hits if "pred" not in x.lower() and "cofold" not in x.lower()]
    _CRYSTAL_CACHE[h] = hits[0] if hits else None
    return _CRYSTAL_CACHE[h]


def main():
    models = [m for m in os.listdir(COFOLD) if os.path.isdir(os.path.join(COFOLD, m))]
    rows, skips = [], []
    meta = {r.holo_id: r for r in HOLOS.itertuples()}
    for _, hr in HOLOS.iterrows():
        h = hr["holo_id"]; lig_code = hr["lig_code"]
        smi = hr.get("smiles_canonical") or hr.get("smiles")
        exp = expected_lig(smi)
        cpath = crystal_path(h)
        if exp is None or cpath is None:
            skips.append((h, "no smiles/crystal")); continue
        n_heavy, ref_el = exp
        try:
            cr = pl.parse_crystal(cpath, lig_code)
        except Exception as e:
            skips.append((h, f"crystal parse {repr(e)[:40]}")); continue
        if not cr["chains"]:
            skips.append((h, "crystal has no protein chain")); continue
        # pick the chain with the most CA (some crystals split protein across chains)
        chain = max(cr["chains"], key=lambda c: len(cr["chains"][c]))
        cca, cres = cr["chains"][chain], cr["chains_res"][chain]
        # all heavy crystal ligand atoms (rectangular match -> nearest copy)
        cxyz, cel = crystal_lig_atoms(cpath, lig_code)
        # must contain at least one full ligand's worth of each element
        if len(cel) < n_heavy or any(Counter(cel)[e] < c for e, c in ref_el.items()):
            skips.append((h, f"crystal lig atoms insufficient ({len(cel)} heavy vs need {dict(ref_el)})")); continue
        for model in models:
            for pf in sorted(glob.glob(f"{COFOLD}/{model}/{h}_s*.pdb")):
                try:
                    p = pl.parse_pose(pf)
                except Exception:
                    continue
                if len(p["ca"]) < 30 or len(p["lig_elems"]) == 0:
                    continue
                al = pl.align_by_sequence(p["ca"], p["ca_res"], cca, cres)
                if al is None or al[2] > 3.0:
                    skips.append((f"{model}/{h}", "align fail")); continue
                R, t, prms, _ = al
                hx, he = heavy_only(p["lig_xyz"], list(p["lig_elems"]))
                samples = split_copies(hx, he, n_heavy, ref_el)
                if not samples:
                    skips.append((f"{model}/{h}", f"pose lig split fail (heavy={len(he)}, need mult of {n_heavy})")); continue
                bf = np.asarray(p["lig_bfac"]) if p.get("lig_bfac") is not None else None
                # OpenProtein returns the ligand duplicated across the N samples; dedupe
                uniq, seen = [], set()
                for sx, se in samples:
                    key = tuple(np.round(sx.mean(0), 2))
                    if key not in seen:
                        seen.add(key); uniq.append((sx, se))
                samples = uniq
                for si, (sx, se) in enumerate(samples):
                    sxT = pl.apply_T(R, t, sx)
                    rmsd = elem_rmsd(sxT, se, cxyz, cel)
                    plddt = float(np.nanmean(bf[si*n_heavy:(si+1)*n_heavy])) if bf is not None and len(bf) >= (si+1)*n_heavy else np.nan
                    rows.append({"model": model, "holo": h, "sample": si, "rmsd": round(rmsd, 3),
                                 "plddt": round(plddt, 2) if not np.isnan(plddt) else np.nan,
                                 "n_heavy": n_heavy, "prot_rmsd": round(prms, 2)})
    df = pd.DataFrame(rows)
    os.makedirs("data/processed/posthoc", exist_ok=True)
    df.to_csv("data/processed/posthoc/cofold_scores.csv", index=False)
    print(f"\n=== scored {len(df)} poses across {df.holo.nunique() if len(df) else 0} holos, {df.model.nunique() if len(df) else 0} models ===")
    if len(df):
        # per (model,holo): best sample = model's realistic self-selected? report both mean-sample and best-sample
        g = df.groupby(["model", "holo"]).rmsd
        best = g.min().groupby("model")
        mean = df.groupby("model").rmsd
        print("\nper-model TRUE RMSD (A):")
        print(f"{'model':14s} {'n_holo':>6s} {'mean_all':>9s} {'best_sample_median':>18s} {'frac<2A(best)':>13s}")
        for m in sorted(df.model.unique()):
            sub = df[df.model == m]
            bh = sub.groupby("holo").rmsd.min()
            print(f"{m:14s} {sub.holo.nunique():6d} {sub.rmsd.mean():9.2f} {bh.median():18.2f} {(bh<2).mean():13.2f}")
        # cross-model oracle per holo (min over all model best-samples)
        holo_best = df.groupby("holo").rmsd.min()
        print(f"\ncross-model ORACLE (best pose any model/sample per holo): median {holo_best.median():.2f}A, frac<2A {(holo_best<2).mean():.2f}")
        # selection wall test: pick by plddt vs oracle
        pk = df.dropna(subset=["plddt"])
        if len(pk):
            sel = pk.loc[pk.groupby("holo").plddt.idxmax()].set_index("holo").rmsd
            orc = pk.groupby("holo").rmsd.min()
            print(f"plddt-select median {sel.median():.2f}A vs oracle {orc.median():.2f}A  (gap {sel.median()-orc.median():+.2f})")
    print(f"\nskips ({len(skips)}):")
    for s in skips[:25]:
        print("  ", s)
    json.dump({"n_scored": len(df), "skips": skips}, open("data/processed/posthoc/cofold_scores_meta.json", "w"))


if __name__ == "__main__":
    main()
