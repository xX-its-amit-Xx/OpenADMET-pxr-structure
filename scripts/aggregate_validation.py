"""Apply re-ranking methods to validation predictions, score each method.

Reads:
  - data/processed/validation_set/scores.csv (per-prediction RMSD vs crystal)
  - For each prediction: its confidence JSON (Boltz output: confidence_*.json)

For each re-ranker method, applies its selection rule per holo to pick one pose,
then reports mean RMSD, median RMSD, success rate (RMSD<2A), and LDDT-PLI proxy.

Output:
  - data/processed/validation_set/method_comparison.csv
  - prints ranked summary table
"""
from __future__ import annotations
import json
import csv
from pathlib import Path
from collections import defaultdict
import numpy as np

REPO    = Path(__file__).resolve().parent.parent
VAL_DIR = REPO / "data" / "processed" / "validation_set"
SCORES  = VAL_DIR / "scores.csv"
OUT_CSV = VAL_DIR / "method_comparison.csv"


def load_scores():
    """Return {holo_id: [(pred_path, rmsd, lddt_proxy), ...]}."""
    by_holo = defaultdict(list)
    with SCORES.open() as f:
        for row in csv.DictReader(f):
            by_holo[row["holo_id"]].append({
                "pred_path": row["pred_path"],
                "lig_rmsd": float(row["lig_rmsd"]),
                "lddt_pli_proxy": float(row["lddt_pli_proxy"]),
            })
    return by_holo


def attach_confidences(by_holo):
    """For each prediction, find its sibling confidence_*.json and attach signals."""
    for holo_id, preds in by_holo.items():
        for p in preds:
            pdb = REPO / p["pred_path"]
            stem = pdb.stem
            # Boltz confidence file: confidence_<stem>.json in same dir
            conf_path = pdb.parent / f"confidence_{stem}.json"
            if not conf_path.exists():
                p["conf"] = {}
                continue
            try:
                p["conf"] = json.loads(conf_path.read_text())
            except Exception:
                p["conf"] = {}


# Selection rules: each takes list-of-preds for one holo, returns the chosen pred dict.

def rule_ipde_best(preds):
    """#262: pick max -complex_ipde (i.e., min complex_ipde)."""
    scored = [(p, -p["conf"].get("complex_ipde", 1e9)) for p in preds]
    scored = [(p, s) for p, s in scored if s > -1e9]
    if not scored:
        return None
    return max(scored, key=lambda x: x[1])[0]


def rule_iplddt_best(preds):
    """#290: pick max complex_iplddt."""
    scored = [(p, p["conf"].get("complex_iplddt", -1)) for p in preds]
    scored = [(p, s) for p, s in scored if s > -1]
    if not scored:
        return None
    return max(scored, key=lambda x: x[1])[0]


def rule_lig_iptm_best(preds):
    """#13 family: pick max ligand_iptm."""
    scored = [(p, p["conf"].get("ligand_iptm", -1)) for p in preds]
    scored = [(p, s) for p, s in scored if s > -1]
    if not scored:
        return None
    return max(scored, key=lambda x: x[1])[0]


def rule_confidence_best(preds):
    """Boltz's headline confidence_score (replaces older sample_ranking_score)."""
    scored = [(p, p["conf"].get("confidence_score", -1)) for p in preds]
    scored = [(p, s) for p, s in scored if s > -1]
    if not scored:
        return None
    return max(scored, key=lambda x: x[1])[0]


def rule_pde_best(preds):
    """Pick min complex_pde (predicted error)."""
    scored = [(p, -p["conf"].get("complex_pde", 1e9)) for p in preds]
    scored = [(p, s) for p, s in scored if s > -1e9]
    if not scored:
        return None
    return max(scored, key=lambda x: x[1])[0]


def rule_plddt_best(preds):
    """Pick max complex_plddt."""
    scored = [(p, p["conf"].get("complex_plddt", -1)) for p in preds]
    scored = [(p, s) for p, s in scored if s > -1]
    if not scored:
        return None
    return max(scored, key=lambda x: x[1])[0]


def rule_oracle(preds):
    """Lower-bound: always pick the lowest-RMSD pose. Tells us the ceiling each rule could hit."""
    return min(preds, key=lambda p: p["lig_rmsd"])


def rule_first(preds):
    """Naive baseline: pick first prediction (model_0)."""
    return preds[0] if preds else None


def rule_top3_median(preds):
    """#300: among top-3 IPDE-best, pick the one closest to atomic median."""
    if not preds:
        return None
    scored = sorted(preds, key=lambda p: p["conf"].get("complex_ipde", 1e9))[:3]
    return scored[len(scored) // 2]  # median by IPDE rank as proxy for atomic median


METHODS = {
    "ORACLE_lower_bound": rule_oracle,
    "#first_model_0": rule_first,
    "#262_ipde_best": rule_ipde_best,
    "#290_iplddt_best": rule_iplddt_best,
    "#13_lig_iptm_best": rule_lig_iptm_best,
    "#conf_score_best": rule_confidence_best,
    "#pde_best": rule_pde_best,
    "#plddt_best": rule_plddt_best,
    "#300_top3_median_proxy": rule_top3_median,
}


def main():
    if not SCORES.exists():
        print(f"No {SCORES} yet — run scripts/score_validation.py first.")
        return
    by_holo = load_scores()
    print(f"Loaded {sum(len(v) for v in by_holo.values())} preds across {len(by_holo)} holos")
    attach_confidences(by_holo)

    # Sanity: how many holos have any confidence data?
    n_with_conf = sum(1 for v in by_holo.values() if any(p["conf"] for p in v))
    print(f"Holos with confidence data: {n_with_conf}/{len(by_holo)}")

    rows = []
    for name, rule in METHODS.items():
        picks = []
        n_picked = 0
        for holo_id, preds in by_holo.items():
            pick = rule(preds)
            if pick:
                picks.append(pick["lig_rmsd"])
                n_picked += 1
        if not picks:
            print(f"  {name}: 0 picks (probably missing confidence field)")
            continue
        arr = np.array(picks)
        rows.append({
            "method": name,
            "n_holos": n_picked,
            "mean_rmsd": round(float(arr.mean()), 3),
            "median_rmsd": round(float(np.median(arr)), 3),
            "frac_under_2A": round(float((arr < 2.0).mean()), 3),
            "frac_under_4A": round(float((arr < 4.0).mean()), 3),
            "mean_lddt_proxy": round(float((1.0 / (1.0 + arr)).mean()), 4),
        })
    rows.sort(key=lambda r: r["mean_rmsd"])

    print("\nRanked by mean_rmsd (lower is better):")
    print(f"  {'method':<28} n  mean_rmsd  median  <2A   <4A   lddt_proxy")
    for r in rows:
        print(f"  {r['method']:<28} {r['n_holos']:>2}  "
              f"{r['mean_rmsd']:>8.3f}  {r['median_rmsd']:>6.3f}  "
              f"{r['frac_under_2A']:>4.2f}  {r['frac_under_4A']:>4.2f}  "
              f"{r['mean_lddt_proxy']:>7.4f}")

    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
