#!/usr/bin/env python3
"""
Prep data for the 3D pose-refinement animation on holo 8CH8 (ligand ULC).
Uses the AUTHORITATIVE scoring path: parse_crystal + align_by_sequence (residue-identity
offset) to put every Boltz pose ligand into the true crystal frame, element-matched RMSD.
Builds a decomposed refinement trajectory (translate->rotate->morph) with per-frame RMSD
for red/green glow. Emits data/processed/posthoc/anim_8CH8.json.
"""
import os
import re
import sys
import json
import glob
import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, "scripts")
import pose_lib as pl

HOLO = "8CH8"; LIGCODE = "ULC"
CRYSTAL = "data/external/novel_val/crystals/8ch8.pdb"


def numsort(paths):
    return sorted(paths, key=lambda p: int(re.search(r"model_(\d+)", p).group(1)))


def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Q.mean(0) - R @ P.mean(0)


def elem_assign(pxT, pe, cx, ce):
    """element-constrained optimal atom assignment pose->crystal. returns idx into crystal per pose atom."""
    D = np.linalg.norm(pxT[:, None, :] - cx[None, :, :], axis=2)
    big = D.max() + 1000.0
    C = D.copy(); C[np.array(pe)[:, None] != np.array(ce)[None, :]] = big
    r, c = linear_sum_assignment(C)
    order = np.zeros(len(pe), dtype=int)
    order[r] = c
    return order, float(np.sqrt((D[r, c] ** 2).mean()))


def crystal_protein_pdb(path):
    return "\n".join(ln.rstrip("\n") for ln in open(path) if ln[:4] == "ATOM")


def main():
    crystal = pl.parse_crystal(CRYSTAL, LIGCODE)
    chain = next(iter(crystal["chains"]))
    cca, cres = crystal["chains"][chain], crystal["chains_res"][chain]
    cp = crystal["lig_copies"][0]
    cx, ce = cp["xyz"], list(cp["elems"])
    print(f"crystal: chain {chain}, {len(cca)} CA, {len(cx)} {LIGCODE} atoms")

    poses = []
    for pf in numsort(glob.glob(f"data/processed/validation_set/preds_nb16/{HOLO}/*.pdb")):
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
    print(f"{len(poses)} poses aligned into crystal frame")
    rmsds = [po["rmsd"] for po in poses]
    oracle = int(np.argmin(rmsds))
    sel = int(np.argmax([po["plddt"] for po in poses]))
    print(f"oracle #{oracle}={rmsds[oracle]:.2f}A | confidence pick #{sel}={rmsds[sel]:.2f}A")

    # crystal atoms reordered to pose-atom order (use oracle's mapping as the canonical map)
    cmap = cx[poses[oracle]["order"]]  # crystal coords in pose-atom order

    # ---- refinement trajectory: confidence pick -> crystal, decomposed ----
    S = poses[sel]["lig"]           # confidence pick ligand (pose-atom order)
    Ccanon = cx[poses[sel]["order"]]  # crystal in this pose's atom order
    frames, edits = [], []

    def rmsd(A, B): return float(np.sqrt(np.mean(np.sum((A - B) ** 2, 1))))

    def add(name, seq):
        for f in seq:
            frames.append([[round(v, 3) for v in a] for a in f.tolist()])
        edits.append({"name": name, "start_frame": len(frames) - len(seq),
                      "end_frame": len(frames) - 1,
                      "rmsd_start": round(rmsd(seq[0], Ccanon), 2),
                      "rmsd_end": round(rmsd(seq[-1], Ccanon), 2)})

    N = 14
    add("start: confidence pick", [S])
    d = Ccanon.mean(0) - S.mean(0)
    add("translate into pocket", [S + d * (k / N) for k in range(N + 1)])
    P1 = S + d
    R2, _ = kabsch(P1, Ccanon)
    ang = np.arccos(np.clip((np.trace(R2) - 1) / 2, -1, 1))
    ax = np.array([R2[2, 1] - R2[1, 2], R2[0, 2] - R2[2, 0], R2[1, 0] - R2[0, 1]])
    ax = ax / (np.linalg.norm(ax) + 1e-9)
    cen = P1.mean(0)
    seq2 = []
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

    out = {
        "holo": HOLO, "ligcode": LIGCODE,
        "protein_pdb": crystal_protein_pdb(CRYSTAL),
        "crystal_lig": [[round(v, 3) for v in a] for a in cmap.tolist()],
        "elems": poses[oracle]["el"],
        "poses": [{"lig": [[round(v, 3) for v in a] for a in po["lig"].tolist()],
                   "rmsd": round(po["rmsd"], 2), "plddt": round(po["plddt"], 1)} for po in poses],
        "oracle_idx": oracle, "selected_idx": sel,
        "traj_frames": frames, "edits": edits, "worst_atoms": worst_atoms,
        "stats": {"n_poses": len(poses), "spread": round(float(np.std(rmsds)), 2),
                  "oracle_rmsd": round(rmsds[oracle], 2), "selected_rmsd": round(rmsds[sel], 2),
                  "pool_min": round(min(rmsds), 2), "pool_max": round(max(rmsds), 2)},
    }
    os.makedirs("data/processed/posthoc", exist_ok=True)
    json.dump(out, open("data/processed/posthoc/anim_8CH8.json", "w"))
    print(f"wrote anim_8CH8.json: {len(frames)} frames, {len(poses)} poses, {len(cx)} atoms | "
          f"oracle {rmsds[oracle]:.2f} sel {rmsds[sel]:.2f}")


if __name__ == "__main__":
    main()
