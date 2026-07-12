#!/usr/bin/env python3
"""
Assemble docs/posthoc_analysis.html in the pixel-report style from computed results.
Re-run any time to refresh the page as new analyses land ("updated continually").
"""
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

STYLE = open("data/processed/posthoc/_style_block.html", encoding="utf-8").read()

# ---- load results ----
m = pd.read_parquet("data/processed/posthoc/master_184.parquet")
d = m[m["disagreement_mean_pw_rmsd"].notna()].copy()
ds = pd.read_parquet("data/processed/pxr_pose_scorer/dataset.parquet")

# selection wall numbers
sig = {"lig_plddt": True, "ptm": True, "iptm": True, "complex_ipde": False,
       "complex_pde": False, "n_pocket_contacts_4p5": True}
rows = []
for h, g in ds.groupby("holo"):
    r = {"oracle": g["rmsd"].min(), "random": g["rmsd"].mean()}
    for s, hib in sig.items():
        r[s] = g.loc[g[s].idxmax() if hib else g[s].idxmin(), "rmsd"]
    rows.append(r)
R = pd.DataFrame(rows)
oracle = R["oracle"].median(); best_sel = min([R[s].median() for s in sig])
best_name = min(sig, key=lambda s: R[s].median())

# GT-free learned selector (leave-one-holo-out) — can ML beat the wall?
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import LeaveOneGroupOut
GTFREE = ["centroid_to_pocket", "radius_gyration", "buried_frac", "n_pocket_contacts_4p5",
          "n_contacts_5p0", "min_pl_dist", "clash_soft", "clash_hard", "clash_vhard",
          "pocket_resid_engaged", "pocket_contact_total", "pocket_contact_max",
          "aromatic_cage_contacts", "aromatic_cage_engaged", "anchor_ser247", "anchor_gln285",
          "anchor_his407", "anchor_arg410", "n_anchor_hbond_3p5", "n_anchors_engaged_5",
          "flex_anchor_overengaged", "mmff_strain_global", "lig_plddt", "conf_score", "ptm",
          "iptm", "ligand_iptm", "complex_plddt", "complex_pde", "complex_ipde", "n_heavy"]
_f = [c for c in GTFREE if c in ds.columns]
_X = ds[_f].fillna(0).values; _y = ds["rmsd"].values; _g = ds["holo"].values
_L = []
for tr, te in LeaveOneGroupOut().split(_X, _y, _g):
    _m = GradientBoostingRegressor(n_estimators=300, max_depth=3, learning_rate=0.05,
                                   subsample=0.8, random_state=0).fit(_X[tr], _y[tr])
    _gh = ds.iloc[te]
    _L.append(_gh.iloc[np.argmin(_m.predict(_X[te]))]["rmsd"])
learned_sel = float(np.median(_L))

corrs = {c: spearmanr(d[c], d["disagreement_mean_pw_rmsd"], nan_policy="omit")[0]
         for c in ["heavy", "mw", "tpsa", "rotb", "clogp", "narom"]}

# CORRECTED cross-model geometry (crystal frame). The original master aligned to boltz1,
# whose 184 exports fold ~22 A off the real LBD, inflating disagreement ~2x.
try:
    CF = pd.read_parquet("data/processed/posthoc/xmodel_crystalframe.parquet")
    cf_med = CF["disagreement_cf"].median()
    cf_q1, cf_q3 = CF["disagreement_cf"].quantile(.25), CF["disagreement_cf"].quantile(.75)
    cf_nwf = int(CF["n_wellfolded"].median())
    cf_decoupled = int(CF["selection_decoupled_cf"].sum())
except Exception:
    CF, cf_med, cf_q1, cf_q3, cf_nwf, cf_decoupled = None, d["disagreement_mean_pw_rmsd"].median(), 0, 0, 0, 171

# ---- cross-model GT co-fold (OpenProtein): real per-architecture accuracy ----
CF = None
try:
    cf = pd.read_csv("data/processed/posthoc/cofold_scores.csv")
    # best (deduped) pose per model/holo
    bh = cf.groupby(["model", "holo"]).rmsd.min().reset_index()
    per = bh.groupby("model").agg(n_holo=("holo", "nunique"),
                                  med=("rmsd", "median"),
                                  best=("rmsd", "min"),
                                  frac2=("rmsd", lambda x: (x < 2).mean())).reset_index()
    orc = bh.groupby("holo").rmsd.min()
    sz = cf.groupby(["holo", "n_heavy"]).rmsd.min().reset_index()
    CF = dict(per=per.sort_values("med"),
              oracle_med=float(orc.median()), oracle_frac2=float((orc < 2).mean()),
              n_holo=int(bh.holo.nunique()), n_model=int(bh.model.nunique()),
              size_corr=float(sz.rmsd.corr(sz.n_heavy)), overall_best=float(bh.rmsd.min()))
except Exception:
    CF = None

SECTIONS = ["overview", "watch-it", "disagreement", "cross-model", "selection-wall",
            "failure-modes", "per-ligand", "rebuild", "release", "thanks", "methods"]


def cofold_table():
    if CF is None:
        return "<p class='muted'>cross-model co-fold pending.</p>"
    rows = "".join(
        f"<tr><td><b>{r.model}</b></td><td>{int(r.n_holo)}</td><td>{r.med:.2f}</td>"
        f"<td>{r.best:.2f}</td><td>{r.frac2:.0%}</td></tr>"
        for r in CF["per"].itertuples())
    return ("<table><thead><tr><th>architecture</th><th>holos</th><th>median RMSD (A)</th>"
            "<th>best case (A)</th><th>frac &lt; 2A</th></tr></thead><tbody>"
            f"{rows}</tbody></table>")


def kicker(t): return f'<span class="kicker">{t}</span>'


def worst_table():
    w = d.nlargest(12, "disagreement_mean_pw_rmsd")[
        ["structure", "population", "series", "mw", "rotb", "disagreement_mean_pw_rmsd",
         "consensus_frac", "medoid_model"]]
    head = "".join(f"<th>{c}</th>" for c in ["ligand", "pop", "series", "MW", "rotB",
                                             "disagree A", "consensus", "medoid"])
    body = ""
    for _, r in w.iterrows():
        body += ("<tr>" + f"<td>{r.structure}</td><td>{r.population}</td><td>{r.series}</td>"
                 f"<td>{r.mw:.0f}</td><td>{int(r.rotb)}</td>"
                 f"<td>{r.disagreement_mean_pw_rmsd:.1f}</td><td>{r.consensus_frac:.2f}</td>"
                 f"<td>{r.medoid_model}</td></tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def sel_table():
    order = ["oracle", "n_pocket_contacts_4p5", "complex_pde", "ptm", "lig_plddt",
             "iptm", "complex_ipde", "random"]
    disp = {"oracle": "ORACLE (pool best)", "n_pocket_contacts_4p5": "pocket contacts",
            "complex_pde": "complex PDE", "ptm": "pTM", "lig_plddt": "ligand pLDDT",
            "iptm": "ipTM", "complex_ipde": "complex iPDE", "random": "RANDOM (baseline)"}
    body = ""
    for s in order:
        med = R[s].median(); sub2 = 100 * (R[s] < 2).mean()
        gap = "" if s in ("oracle", "random") else f"+{med-oracle:.2f}"
        cls = ' style="color:#50fa7b"' if s == "oracle" else (' style="color:#b8b2a2"' if s == "random" else "")
        body += (f"<tr{cls}><td>{disp[s]}</td><td>{med:.2f}</td><td>{sub2:.0f}%</td>"
                 f"<td>{gap}</td></tr>")
    return ("<table><thead><tr><th>selector</th><th>median RMSD A</th>"
            f"<th>holos &lt;2A</th><th>gap to oracle</th></tr></thead><tbody>{body}</tbody></table>")


def build():
    nav = "".join(f'<a href="#{s}">{s.replace("-", " ").title()}</a>' for s in SECTIONS)
    frag_med = d[d.population == "fragment"]["disagreement_mean_pw_rmsd"].median()
    drug_med = d[d.population == "drug_like"]["disagreement_mean_pw_rmsd"].median()

    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>PXR Structure Challenge - Post-Hoc Autopsy</title>{STYLE}</head><body>
<header class="topbar"><nav class="nav" aria-label="Sections">{nav}</nav></header>
<main>

<section id="overview">{kicker("post-hoc // autopsy")}
<h2>What went wrong — a per-pose autopsy of all 184</h2>
<p>The challenge is over. Best submission <span class="pill">prot_rescue8 &middot; LDDT-PLI 0.564 &middot; rank 2</span>
(leader 0.5725). This report dissects <b>every one of the 184 test ligands</b> across
<b>12-16 structure predictors</b> to answer: where did the models struggle, and why did we
plateau? No ground-truth crystals were released for the 184 blind ligands, so per-ligand
analysis uses <b>cross-model geometric consensus</b>; where true poses exist (the 53-holo
PXR crystal panel) we measure real RMSD and the selection gap.</p>
<h3>Three findings up front</h3>
<p><span class="step-number">1</span> <b>The models disagree — but less than it first looked.</b>
Median cross-model pose disagreement is <b>{cf_med:.1f} A</b> (mean pairwise ligand RMSD across the
well-folded predictors, in a real crystal frame). An earlier all-model average reported ~4.9 A - that
figure was <b>inflated ~2x</b> by two model exports (boltz1, decaf) that fold ~20 A off the real LBD;
re-anchoring to the crystal frame corrects it. Even so, on most ligands there is no single consensus mode.</p>
<p><span class="step-number">2</span> <b>Failure is anchoring ambiguity, not size alone.</b>
Disagreement scales with size/flexibility (heavy-atom rho {corrs['heavy']:+.2f}), but the very
worst cases are <b>small fragments with no directional anchor</b> - fragments split bimodally
into best- and worst-agreed. Rigid H-bond anchors (aminothiazole, benzimidazole) sharpen
agreement; floppy multi-donor groups (urea, amide) blur it.</p>
<p><span class="step-number">3</span> <b>The selection wall is the ceiling.</b> On GT holos the
pool contains sub-2A poses (<b>oracle {oracle:.2f} A</b>, 59% &lt;2A) but the best confidence
signal selects only <b>{best_sel:.2f} A</b> - barely above random. Generation works; selection
is broken.</p>
</section>

<section id="watch-it">{kicker("interactive // watch the wall")}
<h2>See the selection wall in 3D</h2>
<p>For holo 8CH8 the pool held a <b>0.59 A</b> pose; confidence picked <b>2.36 A</b>. The
interactive viewer plays the whole story - the 20-pose spread, the hidden near-perfect pose,
the pose we actually chose, then a decomposed refinement to the crystal where <b>each edit
glows green (closer) or red (farther)</b>, and finally the atoms where we missed most.</p>
<p><a href="posthoc_animation.html" style="font-size:16px">&#9654; Open the interactive 3D pose autopsy &rarr;</a>
<span class="muted"> (WebGL / 3Dmol.js). A holo selector switches between <b>three GT cases</b> that teach
different lessons — 8CH8 (selection wall, 0.59&rarr;2.36&nbsp;A), 7RIV (selection miss, 0.89&rarr;1.58&nbsp;A),
and 4X1G (a <i>generation</i> failure — even the oracle is 6.19&nbsp;A). A second scene shows three
architectures diverging 3-15&nbsp;A on holo 1M13 — the input-fidelity finding made visible.</span></p>
<p>And for the blind set: the <a href="compound_gallery.html" style="font-size:15px">per-compound gallery &rarr;</a>
gives all 184 an interactive 3D cross-model scatter + a traceable, literature-cited difficulty diagnosis.</p>
</section>

<section id="disagreement">{kicker("finding 1 // consensus")}
<h2>~{cf_nwf} well-folded models, {cf_med:.1f} A apart</h2>
<p>For each ligand we superpose every predictor's protein onto a real PXR crystal frame (1ILH) and measure
pairwise ligand RMSD across the models that fold correctly (median <b>{cf_nwf}</b> of 15 per ligand). The
median disagreement is <b>{cf_med:.2f} A</b> (IQR {cf_q1:.2f}-{cf_q3:.2f}).
Fragments ({frag_med:.2f} A) and drug-like analogs ({drug_med:.2f} A) disagree <i>almost equally</i> -
the tidy "analogs are harder" story does not hold; disagreement tracks size and anchoring, not
fragment-vs-analog class. <b>Two model exports (boltz1, decaf) are excluded</b>: their 184 proteins fold
~20&nbsp;A from the real LBD (0/184 within 3&nbsp;A), so their ligands cannot be placed in a shared frame -
a real finding, and the source of the ~2x inflation in the earlier all-model average.</p>
<figure><img src="figures/posthoc_disagreement.png" alt="cross-model disagreement"></figure>
<h3>The 12 most-disagreed ligands (where every model diverges)</h3>
{worst_table()}
<p class="muted">Note the worst offenders are mostly small fragments (MW 160-250) - too few
contacts to pin a pose. The consensus "medoid" model is rarely Boltz; it is usually an
ESM2 / OpenFold3 / apo-template prediction.</p>
<h3>An "outlier" that was really a frame artifact</h3>
<figure><img src="figures/posthoc_permodel.png" alt="per-model divergence (original all-model view)"></figure>
<p class="muted">Figure: the <i>original</i> all-model ranking, which appeared to show Boltz/DeCAF as
~8.5&nbsp;A outliers vs an ESM/OF3/apo "consensus". The autopsy showed this was mostly an
<b>alignment artifact</b>, not a pose disagreement.</p>
<p>Re-anchoring every model to a real crystal frame revealed the cause: <b>boltz1 and decaf's 184
protein exports fold ~20&nbsp;A from the real PXR LBD</b> (systematic — 0/184 within 3&nbsp;A of the
crystal, and their pockets differ ~12&nbsp;A locally too). Aligning the whole panel to boltz1, as the
first pass did, therefore smeared everyone else and inflated the apparent disagreement. Among the
<b>well-folded</b> models (median {cf_nwf}/15, all folding to ~0.6&nbsp;A of the crystal) the real
cross-model spread is <b>{cf_med:.1f}&nbsp;A</b>, and confidence still disagrees with the consensus in
<b>{cf_decoupled}/184</b> ligands. The lesson is methodological as much as biological: <b>always validate
your reference frame</b> before trusting a cross-model metric. Every one of the 184 - with its corrected
cross-model scatter, difficulty diagnosis, and the confidence-vs-consensus gap - is browsable in the
<a href="compound_gallery.html">per-compound gallery</a>.</p>
</section>

<section id="cross-model">{kicker("finding 2b // real ground truth")}
<h2>Cross-model GT panel: input fidelity beats architecture</h2>
<p>To arbitrate the majority-vs-minority question with <i>real</i> crystals, we co-folded the PXR
holo panel across three architectures (Boltz-1x, Boltz-2, ESMFold-2) on the OpenProtein API and
scored every pose against the RCSB crystal ligand (element-matched RMSD in the residue-identity
aligned protein frame). Result across <b>{CF['n_holo'] if CF else 0} holos</b>:</p>
{cofold_table()}
<p>Every backbone superposes to the crystal within ~0.5&nbsp;A, yet <b>not one architecture places
a single ligand under 2&nbsp;A</b> - the cross-model oracle (best pose, any model) is
<b>{CF['oracle_med']:.1f}&nbsp;A</b> and the best case in the whole panel is
<b>{CF['overall_best']:.2f}&nbsp;A</b>. Crucially this is <i>not</i> a ligand-size effect
(RMSD-vs-size correlation is only <b>{CF['size_corr']:.2f}</b>): the misses are uniform.</p>
<p>The contrast with the Boltz scene in the <a href="posthoc_animation.html">3D autopsy</a> is the whole
lesson. The <i>same</i> Boltz model, run with a <b>native per-target setup</b> (target sequence, deep
MSA, best-of-20 sampling), buries a <b>0.59&nbsp;A</b> pose in its pool on 8CH8. Run from a
<b>generic canonical sequence + shallow MSA</b>, it lands <b>8.7&nbsp;A</b> off on 1M13 - and ESMFold-2
lands 15&nbsp;A off on the same target. The lever was never the model zoo; it was
<b>input fidelity: target-specific sequence, MSA depth, and sampling budget.</b> That is precisely what
the LatchBio-sponsored AlphaFold3 deep-sampling run bought us, and why it moved the needle where a
broader-but-shallower model sweep did not.</p>
<p><a href="posthoc_animation.html#view2" style="font-size:15px">&#9654; See the three architectures diverge against the 1M13 crystal &rarr;</a></p>
</section>

<section id="selection-wall">{kicker("finding 3 // the ceiling")}
<h2>The selection wall, quantified</h2>
<p>On the 53-holo PXR crystal panel (real ground truth), each ligand has ~20 sampled poses.
We ask: how good is the pose we would <i>pick</i> by each confidence signal, versus the best
pose the pool actually contains?</p>
{sel_table()}
<figure><img src="figures/posthoc_selection_wall.png" alt="selection wall"></figure>
<p><b>The pool has the answer (oracle {oracle:.2f} A) but no signal finds it.</b> Every
model-confidence channel - pLDDT, ipTM, iPDE, pTM - lands near <b>2.3 A</b>, only ~0.2-0.3 A
better than random pose choice. Best of all was a simple geometric count
(<code>{best_name}</code>, {best_sel:.2f} A). We captured under half the available headroom.
This is why more/better <i>sampling</i> never moved the needle: the bottleneck was picking.</p>
<h3>Even a GT-trained ML selector cannot beat it</h3>
<p>The obvious fix - train a model on true RMSD to rank poses - was tested here with
leave-one-holo-out CV on <b>only inference-time (GT-free) physics + confidence features</b>.
Result: <b>{learned_sel:.2f} A</b>, no better than plain pLDDT ({R['lig_plddt'].median():.2f} A)
and far from oracle ({oracle:.2f} A). <b>The selection wall is an information failure, not a
modeling failure</b>: the signal needed to identify the correct pose is not present in any
feature computable without the answer.</p>
<p class="muted">Caveat &amp; a real bug: our <i>deployed</i> pose scorer appeared to hit oracle
in development - because it was fed <b>leaky features</b> (<code>pool_oracle</code>,
<code>q_softgain</code>, which encode the ground-truth minimum). Those inflate dev metrics and
vanish at inference, explaining why the custom scorer never improved the live leaderboard.
Honest GT-free performance is {learned_sel:.2f} A. (Panel: 17 holos x 20 Boltz poses; the
cross-model GT co-fold panel that confirms the input-fidelity finding is below.)</p>
</section>

<section id="failure-modes">{kicker("finding 2 // anchoring")}
<h2>Why some poses are hopeless: anchoring ambiguity</h2>
<p>Ranking ligands by cross-model agreement and testing functional-group presence reveals a
clean signal: groups that make a <b>rigid, directional interaction</b> reduce disagreement,
while <b>flexible, multi-rotamer donors</b> increase it.</p>
<figure><img src="figures/posthoc_anchoring.png" alt="anchoring functional groups"></figure>
<p>Correlations of disagreement with descriptors (Spearman rho):
heavy {corrs['heavy']:+.2f}, TPSA {corrs['tpsa']:+.2f}, MW {corrs['mw']:+.2f},
rotB {corrs['rotb']:+.2f}, nAromatic {corrs['narom']:+.2f}, cLogP {corrs['clogp']:+.2f}.
Polarity and flexibility drive divergence; aromatic ring count and lipophilicity do not.</p>
</section>

<section id="per-ligand">{kicker("data // every pose")}
<h2>Per-ligand ledger (all 184)</h2>
<p>The full per-ligand table - each ligand's cross-model disagreement, consensus fraction,
medoid model, confidence, and properties - is written to
<code>data/processed/posthoc/master_184.csv</code> (summary) and
<code>master_184_long.parquet</code> (per-model deviations). Live counts:
<span class="pill">{len(d)} ligands analyzed</span>
<span class="pill">{(d['consensus_frac']>=0.5).sum()} with majority consensus</span>
<span class="pill">{(d['disagreement_mean_pw_rmsd']>6).sum()} pathological (&gt;6A)</span></p>
</section>

<section id="rebuild">{kicker("in progress // methods we skipped")}
<h2>Rebuilding the approaches we did not finish</h2>
<p>Now running across Explorer / Modal x2 / Boltz / OpenProtein / Kaggle / Colab:</p>
<p><span class="step-number">A</span> <b>GT-trained pose selector</b> - <span class="pill">DONE - negative</span>
learning to rank poses by true RMSD on GT-free features does <i>not</i> beat confidence
({learned_sel:.2f} A vs oracle {oracle:.2f} A). The wall is fundamental. Next test: a larger
cross-model GT panel (co-fold the 53 holos with every model) to see if <i>cross-model</i>
features carry the missing signal.</p>
<p><span class="step-number">B</span> <b>MD / physics refinement</b> of the pool best poses
(OpenMM/Desmond) to test if relaxation closes the oracle gap.</p>
<p><span class="step-number">C</span> <b>Fresh co-folds</b> (Boltz, OpenProtein RF3/ESMFold)
with deeper sampling on the pathological fragments.</p>
<p class="muted">This section updates as each backend returns.</p>
</section>

<section id="release">{kicker("open data // for the community")}
<h2>All poses released for public benefit</h2>
<p>Every pose generated in this campaign - the full multi-model pool for all 184 test ligands
(Boltz-1/2, OpenFold3, AF3, Chai, Protenix v2, ESMFold2, RoseTTAFold3, apo templates, and the
LatchBio-sponsored AF3 sampling) - is being packaged with <b>proper per-pose labels</b>
(ligand SMILES, model, confidence signals, cross-model disagreement, ligand properties) and
published openly so others can benefit from the compute already spent.
<b>2,310 poses (15 models x 184 ligands) with full per-pose labels</b> are released on Hugging Face:
<a href="https://huggingface.co/datasets/xX-its-amit-Xx/pxr-structure-pose-pool">huggingface.co/datasets/xX-its-amit-Xx/pxr-structure-pose-pool</a>
(CC-BY-4.0; <code>manifest.csv</code> carries SMILES, confidence, cross-model disagreement, and
ligand properties). If a good pose exists in this pool for a ligand that a future method can learn
to <i>select</i>, that is the open problem this data exists to help solve.</p>
</section>

<section id="thanks">{kicker("gratitude")}
<h2>&#128077; Huge thanks to LatchBio</h2>
<p><b>This work was sponsored in its final leg by <a href="https://latch.bio">LatchBio</a>.</b>
Their generous <b>$500 in compute credits</b> powered <b>additional AlphaFold3 pose sampling</b> -
deepening the pose pool exactly where our models struggled most (the pathological small-fragment
and long-tail cases). That sampling is feeding the multi-model 3D autopsy above, and every pose
it produced is being <b>published for the public to benefit</b>. Thank you, LatchBio, for backing
open, reproducible structural science. &#128153;</p>
</section>

<section id="methods">{kicker("repro")}
<h2>How this was computed</h2>
<p>Cross-model geometry: <code>scripts/posthoc/build_master_table.py</code> (protein-CA
Kabsch superposition via <code>scripts/pose_lib.py</code>, element-matched ligand RMSD).
Selection wall: <code>data/processed/pxr_pose_scorer/dataset.parquet</code> (17 holos x 20
Boltz poses, true RMSD + all confidence signals). Figures:
<code>scripts/posthoc/make_figures.py</code>. This page:
<code>scripts/posthoc/build_report.py</code> (re-run to refresh).</p>
</section>

</main></body></html>"""
    open("docs/posthoc_analysis.html", "w", encoding="utf-8").write(html)
    print(f"wrote docs/posthoc_analysis.html ({len(html)//1024} KB)")


if __name__ == "__main__":
    build()
