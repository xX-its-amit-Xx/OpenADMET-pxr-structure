#!/usr/bin/env python3
"""Post-hoc figures in the pixel-report dark theme -> docs/figures/posthoc_*.png"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

BG = "#0b0b12"; PANEL = "#171821"; INK = "#f7f1db"; MUTED = "#b8b2a2"
CYAN = "#00e5ff"; PINK = "#ff3d81"; GREEN = "#50fa7b"; YELLOW = "#ffe66d"
ORANGE = "#ff9f1c"; RED = "#ff5a5f"; LINE = "#3f4256"
plt.rcParams.update({
    "figure.facecolor": BG, "axes.facecolor": PANEL, "savefig.facecolor": BG,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.edgecolor": LINE, "font.family": "monospace", "font.size": 11,
    "axes.grid": True, "grid.color": "#23252f", "grid.linewidth": 0.6,
})
OUT = "docs/figures"

# ---------- Fig 1: the selection wall ----------
df = pd.read_parquet("data/processed/pxr_pose_scorer/dataset.parquet")
signals = {"lig_plddt": True, "conf_score": True, "ptm": True, "iptm": True, "ligand_iptm": True,
           "complex_plddt": True, "complex_ipde": False, "complex_pde": False,
           "mmff_strain_global": False, "clash_hard": False, "buried_frac": True,
           "n_pocket_contacts_4p5": True, "n_anchors_engaged_5": True, "min_pl_dist": False}
res = []
for holo, g in df.groupby("holo"):
    row = {"oracle": g["rmsd"].min(), "random": g["rmsd"].mean()}
    for s, hib in signals.items():
        idx = g[s].idxmax() if hib else g[s].idxmin()
        row[s] = g.loc[idx, "rmsd"]
    res.append(row)
R = pd.DataFrame(res)
labels = ["oracle", "n_pocket_contacts_4p5", "complex_pde", "ptm", "lig_plddt",
          "iptm", "complex_ipde", "random"]
disp = {"oracle": "ORACLE (pool best)", "n_pocket_contacts_4p5": "pocket contacts",
        "complex_pde": "complex PDE", "ptm": "pTM", "lig_plddt": "ligand pLDDT",
        "iptm": "ipTM", "complex_ipde": "complex iPDE", "random": "RANDOM"}
meds = [R[l].median() for l in labels]
colors = [GREEN] + [CYAN] * 6 + [MUTED]
fig, ax = plt.subplots(figsize=(9, 4.6))
y = np.arange(len(labels))[::-1]
ax.barh(y, meds, color=colors, edgecolor=LINE, height=0.62)
ax.axvline(meds[0], color=GREEN, ls="--", lw=1, alpha=0.7)
ax.axvline(2.0, color=RED, ls=":", lw=1.2, alpha=0.8)
for yi, m in zip(y, meds):
    ax.text(m + 0.03, yi, f"{m:.2f}A", va="center", color=INK, fontsize=10)
ax.set_yticks(y); ax.set_yticklabels([disp[l] for l in labels], fontsize=10)
ax.set_xlabel("median ligand RMSD of SELECTED pose  (lower = better)")
ax.set_title("THE SELECTION WALL — pool holds sub-2A poses, no signal finds them",
             color=INK, fontsize=12, loc="left", pad=10)
ax.text(2.02, len(labels) - 0.5, "2A", color=RED, fontsize=9)
ax.set_xlim(0, max(meds) * 1.18)
plt.tight_layout(); plt.savefig(f"{OUT}/posthoc_selection_wall.png", dpi=140); plt.close()
print("wrote posthoc_selection_wall.png")

# ---------- Fig 2: cross-model disagreement ----------
m = pd.read_parquet("data/processed/posthoc/master_184.parquet")
d = m[m["disagreement_mean_pw_rmsd"].notna()]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
frag = d[d["population"] == "fragment"]["disagreement_mean_pw_rmsd"]
drug = d[d["population"] == "drug_like"]["disagreement_mean_pw_rmsd"]
a1.hist([frag, drug], bins=20, stacked=True, color=[CYAN, PINK], edgecolor=BG,
        label=[f"fragment (n={len(frag)})", f"drug-like (n={len(drug)})"])
a1.axvline(d["disagreement_mean_pw_rmsd"].median(), color=YELLOW, ls="--", lw=1.2)
a1.text(d["disagreement_mean_pw_rmsd"].median() + 0.1, a1.get_ylim()[1] * 0.9,
        f"median {d['disagreement_mean_pw_rmsd'].median():.1f}A", color=YELLOW, fontsize=9)
a1.set_xlabel("cross-model disagreement (mean pairwise ligand RMSD, A)")
a1.set_ylabel("ligands"); a1.legend(fontsize=8, facecolor=PANEL, edgecolor=LINE, labelcolor=INK)
a1.set_title("12-16 models disagree ~5A on where each ligand binds", color=INK, fontsize=11, loc="left")
# scatter disagreement vs heavy atoms
a2.scatter(d["heavy"], d["disagreement_mean_pw_rmsd"], c=d["rotb"], cmap="cool", s=26,
           edgecolor=LINE, linewidth=0.3)
a2.set_xlabel("heavy atoms"); a2.set_ylabel("disagreement (A)")
a2.set_title("bigger/flexible -> more disagreement (rho +0.28)", color=INK, fontsize=11, loc="left")
sc = a2.collections[0]; cb = fig.colorbar(sc, ax=a2); cb.set_label("rot. bonds", color=MUTED)
cb.ax.yaxis.set_tick_params(color=MUTED)
plt.tight_layout(); plt.savefig(f"{OUT}/posthoc_disagreement.png", dpi=140); plt.close()
print("wrote posthoc_disagreement.png")

# ---------- Fig 3: anchoring functional groups ----------
fgs = [("fg_aminothiazole", "aminothiazole"), ("fg_benzimidazole", "benzimidazole"),
       ("fg_sulfone", "sulfone"), ("fg_sulfonamide", "sulfonamide"),
       ("fg_amide", "amide"), ("fg_urea", "urea")]
deltas, names, ns = [], [], []
base = d["disagreement_mean_pw_rmsd"].median()
for col, nm in fgs:
    if col in d and d[col].sum() >= 3:
        deltas.append(d[d[col] == 1]["disagreement_mean_pw_rmsd"].median() - d[d[col] == 0]["disagreement_mean_pw_rmsd"].median())
        names.append(nm); ns.append(int(d[col].sum()))
fig, ax = plt.subplots(figsize=(8, 4.2))
cols = [GREEN if x < 0 else ORANGE for x in deltas]
yb = np.arange(len(names))
ax.barh(yb, deltas, color=cols, edgecolor=LINE, height=0.6)
ax.axvline(0, color=INK, lw=1)
for yi, x, nn in zip(yb, deltas, ns):
    ax.text(x + (0.03 if x >= 0 else -0.03), yi, f"{x:+.2f}A (n={nn})",
            va="center", ha="left" if x >= 0 else "right", color=INK, fontsize=9)
ax.set_yticks(yb); ax.set_yticklabels(names)
ax.set_xlabel("delta disagreement vs ligands lacking the group  (<0 = clearer pose)")
ax.set_title("ANCHORING drives agreement: rigid H-bond anchors help, floppy donors hurt",
             color=INK, fontsize=11.5, loc="left", pad=10)
ax.set_xlim(min(deltas) * 1.5 - 0.2, max(deltas) * 1.5 + 0.3)
plt.tight_layout(); plt.savefig(f"{OUT}/posthoc_anchoring.png", dpi=140); plt.close()
print("wrote posthoc_anchoring.png")

# ---------- Fig 4: per-model divergence from consensus (architecture clustering) ----------
lg = pd.read_parquet("data/processed/posthoc/master_184_long.parquet")
dev = lg.dropna(subset=["rmsd_to_medoid"]).groupby("model")["rmsd_to_medoid"].median().sort_values()
fig, ax = plt.subplots(figsize=(8.5, 4.6))
yb = np.arange(len(dev))[::-1]
# color: ESM/consensus cluster cyan, outliers pink
cols = [PINK if v > 4 else (YELLOW if v > 2.4 else CYAN) for v in dev.values]
ax.barh(yb, dev.values, color=cols, edgecolor=LINE, height=0.66)
for yi, v in zip(yb, dev.values):
    ax.text(v + 0.1, yi, f"{v:.1f}A", va="center", color=INK, fontsize=9)
ax.set_yticks(yb); ax.set_yticklabels(dev.index, fontsize=9.5)
ax.set_xlabel("median ligand-pose distance from cross-model consensus (A)")
ax.set_title("Models cluster by ARCHITECTURE, not truth - boltz/decaf are outliers",
             color=INK, fontsize=11.5, loc="left", pad=10)
ax.text(0.98, 0.04, "consensus != correctness (no GT for the 184)", transform=ax.transAxes,
        ha="right", color=MUTED, fontsize=8.5, style="italic")
plt.tight_layout(); plt.savefig(f"{OUT}/posthoc_permodel.png", dpi=140); plt.close()
print("wrote posthoc_permodel.png")
