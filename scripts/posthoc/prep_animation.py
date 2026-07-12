#!/usr/bin/env python3
"""
Prep data for the 3D pose-refinement animation, for ANY holo (default 8CH8).

Usage:  python scripts/posthoc/prep_animation.py [HOLO] [LIGCODE]
Uses the AUTHORITATIVE scoring path: parse_crystal + align_by_sequence (residue-identity
offset) to put every Boltz pose ligand into the true crystal frame, element-matched RMSD.
Builds a decomposed refinement trajectory (translate->rotate->morph) with per-frame RMSD,
and per-holo "where we missed" sites (worst atoms -> nearest crystal pocket residues, with
a data-derived, research-cited explanation). Emits data/processed/posthoc/anim_<HOLO>.json.
"""
import os, re, sys, json, glob
import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, "scripts")
import pose_lib as pl

HOLO = (sys.argv[1] if len(sys.argv) > 1 else "8CH8").upper()
LIGCODE = sys.argv[2] if len(sys.argv) > 2 else {"8CH8": "ULC", "7RIV": "P06", "3HVL": "SRL", "4X1G": "3WF"}.get(HOLO)
CRYSTAL_ROOTS = ["data/external/novel_val/crystals", "data/external/pxr_holo", "data/external/rcsb_holo"]
PREDS = f"data/processed/validation_set/preds_nb16/{HOLO}"

# 8CH8 keeps its deep-research-backed site prose (verbatim from why_we_missed_research.md).
RESEARCHED_SITES_8CH8 = [
  {"name": "Hydrophobic arm", "atoms": [9], "group": "tetrahydronaphthalene ring", "res": "Ile414 · Phe420 · Leu411", "disp": 6.6, "color": "#ffe66d",
   "why": "PXR's pocket is ~1,150 Å³ and <b>essentially uncharged</b> — 20 of ~28 lining residues are hydrophobic — so this bicyclic ring gets <b>no directional contact</b> to orient it and can sit in any of five near-isoenergetic sub-sites. And AF3-class models are biased to pin <i>polar</i> head-groups while <b>actively mis-placing hydrophobic groups</b> — ligand-orientation-within-pocket is their single biggest failure mode (58.6% of errors) — so the ring landed in the wrong sub-site, 6.6 Å off.",
   "cite": "Sci Rep 2018 (uncharged pocket, 5 SR12813 orientations); AF3 docking-failure analysis 2025."},
  {"name": "Polar anchor", "atoms": [10, 5], "group": "thiocarbamate C=S", "res": "Ser247 · Gln285 · Cys284", "disp": 4.2, "color": "#ff5a5f",
   "why": "Ser247/Gln285/His407 are the <b>polar anchor cluster</b> that pins canonical agonists like SR12813, and Cys284 is part of the same hydrogen-bond hot-spot. But our ligand only offers a <b>thiocarbonyl (C=S)</b> here, and sulfur is a <b>weak, anisotropic H-bond acceptor</b> — ~20% weaker than a carbonyl O and accepting only along specific lone-pair directions. That shallow, directional signal is easy to miss, so it slipped ~4 Å.",
   "cite": "Sci Rep 2018 (S247/Q285/H407 anchoring); hot-spot H-bond map; ACS/PMC (C=S weak acceptor)."},
  {"name": "His407 contact", "atoms": [0, 1, 2], "group": "methoxypyridine O", "res": "His407", "disp": 4.5, "color": "#ff9f1c",
   "why": "His407 is a key <b>directional anchor</b> — it accepts a hydrogen bond from agonist hydroxyls (T0901317, SR12813). Here the ligand only reaches it with a <b>weak edge-on methoxy contact</b>, and His407's side chain is mobile in a pocket that reshapes around each ligand, so the model placed the group ~4.5 Å away.",
   "cite": "Cell/Structure 2023 (H407↔hydroxyl H-bond); hot-spot map (Gln285·Ser247·His407·Cys284)."},
]


def numsort(paths):
    return sorted(paths, key=lambda p: int(re.search(r"model_(\d+)", p).group(1)))


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Q.mean(0) - R @ P.mean(0)


def elem_assign(pxT, pe, cx, ce):
    D = np.linalg.norm(pxT[:, None, :] - cx[None, :, :], axis=2)
    big = D.max() + 1000.0
    C = D.copy(); C[np.array(pe)[:, None] != np.array(ce)[None, :]] = big
    r, c = linear_sum_assignment(C)
    order = np.zeros(len(pe), dtype=int); order[r] = c
    return order, float(np.sqrt((D[r, c] ** 2).mean()))


def crystal_path():
    for r in CRYSTAL_ROOTS:
        for name in (HOLO.lower(), HOLO.upper()):
            p = os.path.join(r, f"{name}.pdb")
            if os.path.exists(p):
                return p
    raise SystemExit(f"no crystal for {HOLO}")


def crystal_protein_pdb(path):
    return "\n".join(ln.rstrip("\n") for ln in open(path) if ln[:4] == "ATOM")


def protein_contacts(path, coord):
    """nearest protein atoms (resname+resseq, atomname, dist) to a 3D point."""
    best = []
    for ln in open(path):
        if ln[:4] != "ATOM":
            continue
        p = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        best.append((float(np.linalg.norm(coord - p)), ln[17:20].strip() + ln[22:26].strip(), ln[12:16].strip()))
    best.sort()
    return best[:4]


def generic_sites(sel_lig, crystal_lig, worst_atoms, elems, cpath):
    """Data-derived sites for non-8CH8 holos: group worst atoms by nearest pocket residue,
    explain via the residue's role (polar anchor vs hydrophobic wall), cite general PXR claims."""
    colors = ["#ffe66d", "#ff5a5f", "#ff9f1c", "#a6f56b"]
    # cluster worst atoms by spatial proximity (simple: greedy within 4A)
    pts = np.array([sel_lig[i] for i in worst_atoms])
    used = [False] * len(worst_atoms); clusters = []
    for i in range(len(worst_atoms)):
        if used[i]:
            continue
        grp = [i]; used[i] = True
        for j in range(i + 1, len(worst_atoms)):
            if not used[j] and np.linalg.norm(pts[i] - pts[j]) < 4.5:
                grp.append(j); used[j] = True
        clusters.append(grp)
    sites = []
    for k, grp in enumerate(clusters[:4]):
        atoms = [worst_atoms[g] for g in grp]
        cen = np.mean([crystal_lig[a] for a in atoms], axis=0)
        con = protein_contacts(cpath, cen)
        disp = float(np.mean([np.linalg.norm(np.array(sel_lig[a]) - np.array(crystal_lig[a])) for a in atoms]))
        resnames = []
        for _, rn, _ in con[:3]:
            r3 = rn[:3]
            if r3 not in [x[:3] for x in resnames]:
                resnames.append(rn)
        top = con[0][1]; top3 = top[:3]
        els = "".join(sorted(set(elems[a] for a in atoms)))
        POLAR3 = {"SER", "GLN", "HIS", "CYS", "ARG", "THR", "TYR", "ASN", "ASP", "GLU", "LYS"}
        if top3 in POLAR3:   # classify by the nearest (displayed) residue, consistently
            why = (f"This group sits against a <b>polar residue ({top})</b> in PXR's sparse anchor set "
                   f"(Ser247/Gln285/His407 pin canonical agonists). The contact here is weak/edge-on, and PXR's "
                   f"pocket is a large, flexible, mostly-uncharged cavity — so the model under-weighted it and the "
                   f"atoms landed ~{disp:.1f} Å off.")
            cite = "Sci Rep 2018 (anchor cluster S247/Q285/H407, uncharged pocket)."
            color = "#ff5a5f"
        else:
            why = (f"This group packs into a <b>hydrophobic sub-pocket (near {top})</b> with <b>no directional "
                   f"contact</b> to orient it — PXR offers five near-isoenergetic sub-sites, and AF3-class models "
                   f"actively mis-place hydrophobic moieties. The atoms sit ~{disp:.1f} Å from the crystal.")
            cite = "Sci Rep 2018 (five subsites); AF3 docking-failure analysis 2025 (hydrophobic mis-placement)."
            color = colors[k % len(colors)]
        sites.append({"name": f"Site {k+1} ({top3.title()})", "atoms": atoms, "group": f"{els} atoms",
                      "res": " · ".join(r[:3].title() + r[3:] for r in resnames[:3]), "disp": round(disp, 1),
                      "color": color, "why": why, "cite": cite})
    return sites


def main():
    cpath = crystal_path()
    crystal = pl.parse_crystal(cpath, LIGCODE)
    chain = max(crystal["chains"], key=lambda c: len(crystal["chains"][c]))
    cca, cres = crystal["chains"][chain], crystal["chains_res"][chain]
    cp = crystal["lig_copies"][0]
    cx, ce = np.asarray(cp["xyz"]), list(cp["elems"])
    print(f"{HOLO}/{LIGCODE}: chain {chain}, {len(cca)} CA, {len(cx)} lig atoms")

    poses = []
    for pf in numsort(glob.glob(f"{PREDS}/*.pdb")):
        p = pl.parse_pose(pf)
        al = pl.align_by_sequence(p["ca"], p["ca_res"], cca, cres)
        if al is None or al[2] > 2.5:
            continue
        R, t, prms, _ = al
        ligT = pl.apply_T(R, t, p["lig_xyz"])
        pe = list(p["lig_elems"])
        if sorted(pe) != sorted(ce):
            continue
        order, r = elem_assign(ligT, pe, cx, ce)
        poses.append({"lig": ligT, "el": pe, "order": order, "rmsd": r,
                      "plddt": float(np.nanmean(p["lig_bfac"]))})
    if len(poses) < 3:
        raise SystemExit(f"only {len(poses)} poses matched for {HOLO} (ligand mismatch?)")
    rmsds = [po["rmsd"] for po in poses]
    oracle = int(np.argmin(rmsds))
    sel = int(np.argmax([po["plddt"] for po in poses]))
    print(f"{len(poses)} poses | oracle #{oracle}={rmsds[oracle]:.2f}A | conf pick #{sel}={rmsds[sel]:.2f}A")

    # canonical order = selected pose's atom order (consistent trajectory + sites)
    cmap = cx[poses[sel]["order"]]
    elems = poses[sel]["el"]
    S = poses[sel]["lig"]
    Ccanon = cx[poses[sel]["order"]]
    frames, edits = [], []

    def rmsd(A, B): return float(np.sqrt(np.mean(np.sum((A - B) ** 2, 1))))

    def add(name, seq):
        for f in seq:
            frames.append([[round(v, 3) for v in a] for a in f.tolist()])
        edits.append({"name": name, "start_frame": len(frames) - len(seq), "end_frame": len(frames) - 1,
                      "rmsd_start": round(rmsd(seq[0], Ccanon), 2), "rmsd_end": round(rmsd(seq[-1], Ccanon), 2)})

    N = 14
    add("start: confidence pick", [S])
    d = Ccanon.mean(0) - S.mean(0)
    add("translate into pocket", [S + d * (k / N) for k in range(N + 1)])
    P1 = S + d
    R2, _ = kabsch(P1, Ccanon)
    ang = np.arccos(np.clip((np.trace(R2) - 1) / 2, -1, 1))
    ax = np.array([R2[2, 1] - R2[1, 2], R2[0, 2] - R2[2, 0], R2[1, 0] - R2[0, 1]])
    ax = ax / (np.linalg.norm(ax) + 1e-9)
    cen = P1.mean(0); seq2 = []
    for k in range(N + 1):
        th = ang * k / N
        K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
        Rk = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)
        seq2.append((Rk @ (P1 - cen).T).T + cen)
    add("rotate core into register", seq2)
    P2 = seq2[-1]
    add("fix torsions / local geometry", [P2 + (Ccanon - P2) * (k / N) for k in range(N + 1)])

    per_atom_miss = np.linalg.norm(P2 - Ccanon, axis=1)
    worst_atoms = [int(i) for i in np.argsort(-per_atom_miss)[:6]]

    sites = RESEARCHED_SITES_8CH8 if HOLO == "8CH8" else generic_sites(S, Ccanon, worst_atoms, elems, cpath)

    out = {
        "holo": HOLO, "ligcode": LIGCODE,
        "protein_pdb": crystal_protein_pdb(cpath),
        "crystal_lig": [[round(v, 3) for v in a] for a in cmap.tolist()],
        "elems": elems,
        "poses": [{"lig": [[round(v, 3) for v in a] for a in po["lig"].tolist()],
                   "rmsd": round(po["rmsd"], 2), "plddt": round(po["plddt"], 1)} for po in poses],
        "oracle_idx": oracle, "selected_idx": sel,
        "traj_frames": frames, "edits": edits, "worst_atoms": worst_atoms, "sites": sites,
        "stats": {"n_poses": len(poses), "spread": round(float(np.std(rmsds)), 2),
                  "oracle_rmsd": round(rmsds[oracle], 2), "selected_rmsd": round(rmsds[sel], 2),
                  "pool_min": round(min(rmsds), 2), "pool_max": round(max(rmsds), 2)},
    }
    os.makedirs("data/processed/posthoc", exist_ok=True)
    p = f"data/processed/posthoc/anim_{HOLO}.json"
    json.dump(out, open(p, "w"))
    print(f"wrote {p}: {len(frames)} frames, {len(poses)} poses, {len(sites)} sites | "
          f"oracle {rmsds[oracle]:.2f} sel {rmsds[sel]:.2f}")


if __name__ == "__main__":
    main()
