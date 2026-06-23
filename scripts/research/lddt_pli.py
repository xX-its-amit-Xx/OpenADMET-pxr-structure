"""
LDDT-PLI (protein-ligand interface lDDT) scorer for the GT harness.

WHY: the leaderboard scores LDDT-PLI (a contact / distance-difference metric that
tolerates rigid-body shifts of the whole ligand as long as the local protein<->ligand
contact geometry is preserved), but our validation harness has been using ligand RMSD.
RMSD can declare a pose "failed" (e.g. 3-4 A) even when the pharmacophore contacts are
essentially correct -- a FALSE NEGATIVE under the real metric. This module implements the
standard OpenStructure / CASP15-16 LDDT-PLI and re-evaluates prior RMSD-based conclusions.

DEFINITION (implemented):
  Reference = crystal complex. For every (protein-heavy-atom i, ligand-heavy-atom j) pair
  whose distance d_ref(i,j) <= inclusion_radius (default 6 A), the pair is "in the interface".
  In the predicted complex (superposed on the reference by the PROTEIN backbone), compute
  d_pred(i',j') for the CORRESPONDING atoms. The pair is "preserved" at tolerance T if
  |d_pred - d_ref| <= T. LDDT-PLI = mean over T in {0.5,1,2,4} A of the fraction of preserved
  interface pairs. Range 0..1, higher better. It is a distance-DIFFERENCE metric: after the
  protein superposition, NO further ligand alignment is done (that is the whole point -- it
  measures whether contacts are right, not whether the ligand is rigid-body placed right).

CORRESPONDENCE:
  - protein atoms: matched by (aligned residue id, atom name).
  - ligand atoms: matched predicted<->crystal by the element-constrained assignment already
    used in pose_lib.gt_min_rmsd (Hungarian on the protein-frame distance matrix), reused here.

Usage:
  python scripts/research/lddt_pli.py sanity   # crystal-vs-self ~1.0 sanity check
  python scripts/research/lddt_pli.py score     # score all 18-holo poses, write json + report
"""
from __future__ import annotations
import os, sys, json, glob, tempfile
from collections import defaultdict
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import pose_lib as PL
from pose_lib import AA3, align_by_sequence, apply_T
from boltz_gt_validate import HC, _cif_to_pdb

EXP = os.path.join(REPO, "data", "external", "boltz_api", "boltz-experiments")
CRYS = os.path.join(REPO, "data", "external", "val_crystals")
GT_SCORES = os.path.join(REPO, "data", "external", "boltz_api", "gt_scores.json")
OUT = os.path.join(REPO, "data", "processed", "lddt_pli_eval.json")

THRESHOLDS = (0.5, 1.0, 2.0, 4.0)


# ----------------------------------------------------------------- crystal heavy atoms
def parse_crystal_full(path, lig_code):
    """Crystal parse keeping ALL protein heavy atoms (LDDT-PLI needs them), CA for
    alignment, and ligand copies. Per chain so homodimers are handled.

    Returns dict:
      chains:      {chain -> {resseq -> CA xyz}}   (for align_by_sequence)
      chains_res:  {chain -> {resseq -> 1letter}}
      prot_atoms:  {chain -> list of (resseq, atom_name, element, xyz)}  heavy protein atoms
      lig_copies:  list of dict(elems, xyz)  (one per ligand instance, any chain)
    """
    chains = defaultdict(dict); chains_res = defaultdict(dict)
    prot_atoms = defaultdict(list); ligs = defaultdict(list)
    one = PL.AA1
    with open(path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            resn = line[17:20].strip(); aname = line[12:16].strip(); chain = line[21]
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            elem = (line[76:78].strip() or aname[:1]).upper()
            if rec == "ATOM" and resn in AA3:
                if elem == "H":
                    continue
                k = int(line[22:26])
                prot_atoms[chain].append((k, aname, elem, np.array([x, y, z])))
                if aname == "CA":
                    chains[chain][k] = np.array([x, y, z])
                    chains_res[chain][k] = one.get(resn, "X")
            elif rec == "HETATM" and resn == lig_code and elem != "H":
                ligs[(chain, line[22:27])].append((elem, [x, y, z]))
    lig_copies = []
    for key in sorted(ligs):
        els = [e for e, _ in ligs[key]]; xs = [c for _, c in ligs[key]]
        lig_copies.append(dict(elems=els, xyz=np.array(xs)))
    return dict(chains=dict(chains), chains_res=dict(chains_res),
                prot_atoms=dict(prot_atoms), lig_copies=lig_copies)


# ----------------------------------------------------------------- core scorer
def _ligand_assignment(pe, pxT, ce, cx):
    """element-constrained Hungarian assignment of predicted ligand atoms (pe,pxT in
    REF frame) to crystal ligand atoms (ce,cx). Returns array map[i]=crystal_index, or None."""
    from scipy.optimize import linear_sum_assignment
    if len(ce) != len(pe) or sorted(ce) != sorted(pe):
        return None
    D = np.linalg.norm(pxT[:, None, :] - cx[None, :, :], axis=2)
    big = D.max() + 1000.0
    mask = (np.array(pe)[:, None] != np.array(ce)[None, :])
    C = D.copy(); C[mask] = big
    r, c = linear_sum_assignment(C)
    out = np.empty(len(pe), dtype=int)
    out[r] = c
    return out


def lddt_pli(pose, crystal, inclusion_radius=6.0, thresholds=THRESHOLDS, return_detail=False):
    """Standard LDDT-PLI for one predicted pose against a crystal reference.

    pose:    pose_lib.parse_pose dict (predicted complex, chain A protein + LIG ligand)
    crystal: parse_crystal_full dict
    Returns float in [0,1] (best over crystal chains / ligand copies), or None if
    no valid superposition / ligand-element match. If return_detail, returns
    (score, dict(n_pairs, prot_rms, chain, rmsd_to_that_copy)).
    """
    pe = np.array(pose["lig_elems"]); px = pose["lig_xyz"]
    if px.shape[0] == 0 or not crystal["lig_copies"]:
        return (None, {}) if return_detail else None
    best = None
    for chain, cca in crystal["chains"].items():
        al = align_by_sequence(pose["ca"], pose["ca_res"], cca, crystal["chains_res"][chain])
        if al is None:
            continue
        R, t, prms, _ = al
        if prms > 2.5:
            continue
        pxT = apply_T(R, t, px)  # predicted ligand in REFERENCE (crystal) frame

        # map aligned protein residues: pose resseq -> crystal resseq (the offset that won)
        # rebuild offset from align_by_sequence's matched pairs by re-deriving best offset.
        off = _best_offset(pose["ca"], pose["ca_res"], cca, crystal["chains_res"][chain])
        if off is None:
            continue
        # predicted protein heavy atoms in ref frame, keyed by (crystal_resseq, atom_name)
        pred_prot = _pose_prot_atoms_in_ref(pose, R, t, off)

        # crystal protein heavy atoms keyed by (resseq, atom_name)
        crys_prot = {}
        for (k, an, el, xyz) in crystal["prot_atoms"][chain]:
            crys_prot[(k, an)] = xyz

        for cp in crystal["lig_copies"]:
            ce = cp["elems"]; cx = cp["xyz"]
            amap = _ligand_assignment(pe, pxT, ce, cx)
            if amap is None:
                continue
            # build interface pairs from the CRYSTAL: (crystal protein atom, crystal lig atom)
            # within inclusion radius, then check preservation against predicted atoms.
            score, npairs = _score_pairs(crys_prot, pred_prot, cx, pxT, ce, pe, amap,
                                         inclusion_radius, thresholds)
            if score is None:
                continue
            # also record an rmsd to this copy for monotonicity diagnostics
            rmsd = float(np.sqrt(((pxT - cx[amap]) ** 2).sum(1).mean()))
            cand = (score, dict(n_pairs=npairs, prot_rms=float(prms), chain=chain, rmsd=rmsd))
            if best is None or cand[0] > best[0]:
                best = cand
    if best is None:
        return (None, {}) if return_detail else None
    return best if return_detail else best[0]


def _score_pairs(crys_prot, pred_prot, cx, pxT, ce, pe, amap, radius, thresholds):
    """crys_prot/pred_prot: {(resseq,atom_name)->xyz}. cx crystal lig xyz, pxT pred lig
    xyz (ref frame), amap[i]=crystal lig index for predicted lig atom i."""
    # invert: for predicted lig atom i, crystal counterpart is amap[i].
    # iterate over crystal interface pairs (prot atom key, crystal lig atom index j).
    # need predicted counterparts: protein atom by same key, lig atom = i where amap[i]==j.
    inv = {int(amap[i]): i for i in range(len(amap))}  # crystal lig idx -> predicted lig idx
    shared_keys = [k for k in crys_prot if k in pred_prot]
    if not shared_keys:
        return None, 0
    cp_xyz = np.array([crys_prot[k] for k in shared_keys])
    pp_xyz = np.array([pred_prot[k] for k in shared_keys])
    nthr = len(thresholds)
    total = 0
    preserved = np.zeros(nthr)
    # distances crystal protein-atom <-> crystal ligand-atom
    Dref = np.linalg.norm(cp_xyz[:, None, :] - cx[None, :, :], axis=2)  # (P, Lc)
    Dpred_full = None
    for jp, jc in enumerate(range(cx.shape[0])):
        ipred = inv.get(jc)
        if ipred is None:
            continue
        ref_col = Dref[:, jc]
        in_iface = ref_col <= radius
        if not in_iface.any():
            continue
        # predicted distances for this lig atom: pred protein atoms <-> pred lig atom ipred
        pred_col = np.linalg.norm(pp_xyz - pxT[ipred][None, :], axis=1)
        dd = np.abs(pred_col[in_iface] - ref_col[in_iface])
        total += int(in_iface.sum())
        for ti, T in enumerate(thresholds):
            preserved[ti] += int((dd <= T).sum())
    if total == 0:
        return None, 0
    per_thr = preserved / total
    return float(per_thr.mean()), total


def _best_offset(mob_ca, mob_res, ref_ca, ref_res):
    """re-derive the integer residue-numbering offset that align_by_sequence picks
    (max residue-identity overlap). Returns offset or None."""
    if not mob_ca or not ref_ca:
        return None
    mob_keys = sorted(mob_ca); ref_keys = set(ref_ca)
    best = None
    lo = min(ref_ca) - max(mob_ca); hi = max(ref_ca) - min(mob_ca)
    for off in range(lo, hi + 1):
        n = sum(1 for k in mob_keys
                if (k + off) in ref_keys and mob_res.get(k) == ref_res.get(k + off))
        if best is None or n > best[1]:
            best = (off, n)
    return best[0] if best and best[1] >= 20 else None


def _pose_prot_atoms_in_ref(pose, R, t, off):
    """Transform predicted protein heavy atoms into ref frame, key by (crystal_resseq,
    atom_name). Needs the full protein atom list -> reparse from the pose block is not
    stored, so we recover from the pose PDB path. Instead we stored only CA in parse_pose;
    so we re-read the pose's source. To stay self-contained we attach prot_atoms in
    parse_pose_full below and pass through pose['prot_atoms']."""
    out = {}
    for (k, an, el, xyz) in pose.get("prot_atoms", []):
        out[(k + off, an)] = apply_T(R, t, xyz[None, :])[0]
    return out


def parse_pose_full(path):
    """Like pose_lib.parse_pose but also keeps protein heavy atoms (resseq, name, elem, xyz)
    needed for LDDT-PLI. Reuses parse_pose for ca/ligand fields."""
    base = PL.parse_pose(path)
    prot = []
    with open(path) as f:
        for line in f:
            if line[:6].strip() != "ATOM":
                continue
            resn = line[17:20].strip()
            if resn not in AA3:
                continue
            aname = line[12:16].strip()
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            elem = (line[76:78].strip() or aname[:1]).upper()
            if elem == "H":
                continue
            prot.append((int(line[22:26]), aname, elem, np.array([x, y, z])))
    base["prot_atoms"] = prot
    return base


# ----------------------------------------------------------------- helpers / loaders
def _cif_to_pose(cif_path):
    tf = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
    tf.write(_cif_to_pdb(open(cif_path).read())); tf.close()
    pose = parse_pose_full(tf.name); os.unlink(tf.name)
    return pose


def _sample_metrics(pdir):
    """{sample_k -> {ipde, iptm, plddt}} from metrics.json."""
    out = {}
    mp = os.path.join(pdir, "metrics.json")
    if not os.path.exists(mp):
        return out
    M = json.load(open(mp))
    for k, s in enumerate(M.get("all_sample_results", [])):
        m = s.get("metrics", {})
        out[k] = dict(ipde=m.get("complex_ipde", 9.9),
                      iptm=m.get("iptm", 0.0),
                      plddt=m.get("complex_plddt", 0.0))
    return out


H2_FLAGGED = ["3HVL", "4J5X", "4NY9", "4X1G", "7AXJ", "7N2A", "8R81", "9FZI"]


# ----------------------------------------------------------------- sanity
def sanity():
    print("=== SANITY: crystal-vs-self should be ~1.0 ===")
    for holo in ["6HJ2", "3HVL", "9FZI"]:
        code = HC[holo]
        cp = os.path.join(CRYS, f"{holo}.pdb")
        crystal = parse_crystal_full(cp, code)
        if not crystal["lig_copies"]:
            print(f"{holo}: no crystal ligand"); continue
        # build a "pose" from the crystal itself: protein chain A CA/atoms + ligand copy 0
        pose = _crystal_as_pose(crystal)
        s = lddt_pli(pose, crystal, inclusion_radius=6.0)
        s4 = lddt_pli(pose, crystal, inclusion_radius=4.0)
        print(f"{holo}: self LDDT-PLI(6A)={s}  (4A)={s4}")


def _crystal_as_pose(crystal):
    """Make a pose dict from the crystal's first chain protein + first ligand copy, so we
    can score the crystal against itself (must give ~1.0)."""
    chain = sorted(crystal["chains"])[0]
    ca = {}; ca_res = {}
    for k, xyz in crystal["chains"][chain].items():
        ca[k] = xyz; ca_res[k] = crystal["chains_res"][chain][k]
    prot = [(k, an, el, xyz) for (k, an, el, xyz) in crystal["prot_atoms"][chain]]
    cp = crystal["lig_copies"][0]
    return dict(ca=ca, ca_res=ca_res, prot_atoms=prot,
                lig_elems=cp["elems"], lig_xyz=cp["xyz"].copy(),
                lig_names=["X"] * len(cp["elems"]))


# ----------------------------------------------------------------- full scoring
def score():
    rmsd_ref = {}
    if os.path.exists(GT_SCORES):
        for r in json.load(open(GT_SCORES)):
            rmsd_ref[(r["holo"], r["arm"])] = r

    all_poses = []   # flat list: dict per pose
    per_holo = {}    # holo -> aggregation
    holos = sorted(set(HC) & set(h[:-4] for h in os.listdir(CRYS) if h.endswith(".pdb")))

    for holo in holos:
        code = HC[holo]
        crystal = parse_crystal_full(os.path.join(CRYS, f"{holo}.pdb"), code)
        if not crystal["lig_copies"]:
            print("no crystal ligand", holo, code); continue
        holo_poses = []
        for arm in ("free", "restr"):
            slug = f"gt-{holo}-{arm}-v1"
            pdir = os.path.join(EXP, slug, "outputs", "files", "prediction")
            cifs = sorted(glob.glob(os.path.join(pdir, "sample_*_predicted_structure.cif")))
            if not cifs:
                continue
            mets = _sample_metrics(pdir)
            for cif in cifs:
                k = int(os.path.basename(cif).split("_")[1])
                pose = _cif_to_pose(cif)
                lp6, det = lddt_pli(pose, crystal, inclusion_radius=6.0, return_detail=True)
                lp4 = lddt_pli(pose, crystal, inclusion_radius=4.0)
                rm, n, prms = PL.gt_min_rmsd(pose, crystal)
                m = mets.get(k, {})
                rec = dict(holo=holo, arm=arm, sample=k,
                           lddt_pli_6=lp6, lddt_pli_4=lp4, rmsd=rm,
                           n_pairs=det.get("n_pairs"),
                           ipde=m.get("ipde"), iptm=m.get("iptm"), plddt=m.get("plddt"))
                holo_poses.append(rec); all_poses.append(rec)
        per_holo[holo] = holo_poses
        if holo_poses:
            valid = [p for p in holo_poses if p["lddt_pli_6"] is not None]
            if valid:
                best = max(valid, key=lambda p: p["lddt_pli_6"])
                sel = min(valid, key=lambda p: (p["ipde"] if p["ipde"] is not None else 9.9))
                print(f"{holo}: best-LDDTPLI={best['lddt_pli_6']:.3f} (rmsd {best['rmsd']:.2f}A)  "
                      f"ipde-sel={sel['lddt_pli_6']:.3f} (rmsd {sel['rmsd']:.2f}A)  "
                      f"n={len(valid)}")

    _analyze(all_poses, per_holo, rmsd_ref)


def _spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    rx = _rankdata(x); ry = _rankdata(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def _rankdata(a):
    a = np.asarray(a, float)
    order = a.argsort(kind="mergesort")
    ranks = np.empty(len(a), float)
    sa = a[order]
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and sa[j + 1] == sa[i]:
            j += 1
        avg = (i + j) / 2.0
        ranks[order[i:j + 1]] = avg
        i = j + 1
    return ranks


def _analyze(all_poses, per_holo, rmsd_ref):
    valid = [p for p in all_poses if p["lddt_pli_6"] is not None and p["rmsd"] is not None]
    lp = np.array([p["lddt_pli_6"] for p in valid])
    rm = np.array([p["rmsd"] for p in valid])

    report = {}
    report["n_poses_scored"] = len(valid)

    # (a) RMSD vs LDDT-PLI correlation
    rho_rmsd = _spearman(rm, lp)
    pear = float(np.corrcoef(rm, lp)[0, 1])
    report["rmsd_vs_lddtpli"] = dict(spearman=rho_rmsd, pearson=pear,
                                     note="negative expected (low rmsd <-> high lddt-pli)")
    print("\n=== (a) RMSD vs LDDT-PLI (all poses) ===")
    print(f"  n={len(valid)}  spearman={rho_rmsd:+.3f}  pearson={pear:+.3f}")
    # tolerance: what LDDT-PLI do poses get at various RMSD bins?
    bins = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 6), (6, 99)]
    print("  RMSD bin   -> mean LDDT-PLI (n)")
    tol_rows = []
    for a, b in bins:
        sel = lp[(rm >= a) & (rm < b)]
        if len(sel):
            print(f"   [{a},{b})A   -> {sel.mean():.3f}  (n={len(sel)})")
            tol_rows.append(dict(rmsd_lo=a, rmsd_hi=b, mean_lddtpli=float(sel.mean()), n=int(len(sel))))
    report["tolerance_bins"] = tol_rows

    # (b) H2-flagged "no sub-2A pose" holos -- what LDDT-PLI do their BEST poses get?
    print("\n=== (b) H2-flagged 'NO sub-2A pose' holos under LDDT-PLI ===")
    flagged = []
    for holo in H2_FLAGGED:
        ps = [p for p in per_holo.get(holo, []) if p["lddt_pli_6"] is not None]
        if not ps:
            continue
        best_lp = max(ps, key=lambda p: p["lddt_pli_6"])
        best_rm = min(ps, key=lambda p: p["rmsd"])
        flagged.append(dict(holo=holo,
                            best_lddtpli=best_lp["lddt_pli_6"], best_lddtpli_rmsd=best_lp["rmsd"],
                            best_rmsd=best_rm["rmsd"], best_rmsd_lddtpli=best_rm["lddt_pli_6"]))
        print(f"  {holo}: best-LDDT-PLI={best_lp['lddt_pli_6']:.3f} (its rmsd {best_lp['rmsd']:.2f}A) | "
              f"min-rmsd={best_rm['rmsd']:.2f}A (its lddt-pli {best_rm['lddt_pli_6']:.3f})")
    report["h2_flagged"] = flagged
    if flagged:
        mlp = np.mean([f["best_lddtpli"] for f in flagged])
        print(f"  --> mean best-LDDT-PLI on the 8 'failed' holos = {mlp:.3f}")
        report["h2_flagged_mean_best_lddtpli"] = float(mlp)

    # (c) IPDE-selected vs best-LDDT-PLI (real-metric ranking loss)
    print("\n=== (c) selection: IPDE-selected vs best-LDDT-PLI (real metric) ===")
    sel_rows = []
    for holo, ps in per_holo.items():
        vs = [p for p in ps if p["lddt_pli_6"] is not None]
        if not vs:
            continue
        oracle = max(vs, key=lambda p: p["lddt_pli_6"])
        ipde_sel = min(vs, key=lambda p: (p["ipde"] if p["ipde"] is not None else 9.9))
        sel_rows.append(dict(holo=holo, oracle_lddtpli=oracle["lddt_pli_6"],
                             ipde_sel_lddtpli=ipde_sel["lddt_pli_6"],
                             loss=oracle["lddt_pli_6"] - ipde_sel["lddt_pli_6"]))
    if sel_rows:
        orc = np.mean([r["oracle_lddtpli"] for r in sel_rows])
        ips = np.mean([r["ipde_sel_lddtpli"] for r in sel_rows])
        print(f"  oracle mean LDDT-PLI      = {orc:.3f}")
        print(f"  IPDE-selected mean        = {ips:.3f}")
        print(f"  ranking loss (oracle-ipde)= {orc - ips:.3f}")
        report["selection"] = dict(oracle_mean=float(orc), ipde_sel_mean=float(ips),
                                   ranking_loss=float(orc - ips), per_holo=sel_rows)

    # (d) which cheap signal correlates best with LDDT-PLI (vs with RMSD)?
    print("\n=== (d) signal correlations: vs LDDT-PLI and vs RMSD (within-holo pooled) ===")
    sig_report = {}
    for sig in ("ipde", "iptm", "plddt"):
        xs = np.array([p[sig] for p in valid if p[sig] is not None], float)
        lpv = np.array([p["lddt_pli_6"] for p in valid if p[sig] is not None], float)
        rmv = np.array([p["rmsd"] for p in valid if p[sig] is not None], float)
        if len(xs) < 5:
            continue
        # within-holo z-scored pooling (selection is per-holo)
        wl, wr, ws = _within_holo_pool(valid, sig)
        rho_lp = _spearman(ws, wl); rho_rm = _spearman(ws, wr)
        sig_report[sig] = dict(spearman_vs_lddtpli=rho_lp, spearman_vs_rmsd=rho_rm)
        print(f"  {sig:6s}: vs LDDT-PLI rho={rho_lp:+.3f}   vs RMSD rho={rho_rm:+.3f}")
    report["signal_correlations"] = sig_report

    # (e) oracle / selection reframed vs leader
    if sel_rows:
        report["reframe"] = dict(
            oracle_mean_lddtpli=float(np.mean([r["oracle_lddtpli"] for r in sel_rows])),
            ipde_sel_mean_lddtpli=float(np.mean([r["ipde_sel_lddtpli"] for r in sel_rows])),
            leader=0.5725, prior_rmsd_oracle_estimate=0.5346,
            note="18-holo GT estimate; leaderboard is 184-ligand and includes hard analogs")
        print("\n=== (e) reframe vs leaderboard ===")
        print(f"  18-holo LDDT-PLI oracle    = {report['reframe']['oracle_mean_lddtpli']:.3f}")
        print(f"  18-holo LDDT-PLI ipde-sel  = {report['reframe']['ipde_sel_mean_lddtpli']:.3f}")
        print(f"  leader (184- set)          = 0.5725")

    json.dump(dict(report=report, poses=valid), open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")


def _within_holo_pool(valid, sig):
    """z-score the signal within each holo (selection is per-holo), return pooled
    (lddtpli, rmsd, signal_z) arrays."""
    byho = defaultdict(list)
    for p in valid:
        if p[sig] is not None:
            byho[p["holo"]].append(p)
    L, Rm, S = [], [], []
    for holo, ps in byho.items():
        s = np.array([p[sig] for p in ps], float)
        if s.std() == 0:
            z = s * 0
        else:
            z = (s - s.mean()) / s.std()
        for pi, p in enumerate(ps):
            L.append(p["lddt_pli_6"]); Rm.append(p["rmsd"]); S.append(z[pi])
    return np.array(L), np.array(Rm), np.array(S)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "score"
    {"sanity": sanity, "score": score}.get(mode, score)()
