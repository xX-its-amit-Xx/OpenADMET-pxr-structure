#!/usr/bin/env python3
"""
Traceable per-compound difficulty engine for all 184 PXR structure-track ligands.

Every per-compound statement is either (a) a computed data value from the master
table, or (b) a reference to an adversarially-verified research claim (CLAIMS below,
from the deep-research pass) or a campaign finding (FINDINGS). NOTHING is free-written
per compound — difficulty tags fire from data thresholds and pull cited text. This is
the anti-hallucination backbone: the UI renders only what this engine emits.

Honest scope: the 184 blind ligands have NO crystal ground truth. "Difficulty" here is
diagnosed from ligand descriptors + cross-model behaviour (disagreement, consensus,
confidence-vs-consensus decoupling), each mapped to why it is hard per the literature.
Real crystal RMSD / "where we missed" is only shown for the holo panel (separate script).

Output: data/processed/posthoc/compound_analysis.json  (+ docs copy)
"""
import json, os
import numpy as np
import pandas as pd

M = pd.read_parquet("data/processed/posthoc/master_184.parquet")
ML = pd.read_parquet("data/processed/posthoc/master_184_long.parquet")

# ----------------------------------------------------------------------------- claims
# Curated subset of the 24 adversarially-verified claims (3-vote), verbatim-short.
# id -> {t: short claim text, s: source domain}. Full text in docs/why_we_missed_research.md
CLAIMS = {
 "pocket_size":   {"t": "PXR's pocket is ~1,150–1,600 A^3 and expands to fit ligands from 250 to >800 Da — far larger than other nuclear receptors.", "s": "Sci Rep 2018; Structure 2023; Oxford"},
 "uncharged":     {"t": "The pocket interior is essentially uncharged/hydrophobic — 20 of ~28 lining residues are hydrophobic (E321–R410 salt bridge neutralises charge), so anchoring is weak and directional cues are scarce.", "s": "Sci Rep 2018"},
 "anchor_cluster":{"t": "Only six side chains bind consistently across all PXR structures: polar Ser247·Gln285·His407 and hydrophobic Met243·Trp299·Phe420; SR12813 anchors via S247/Q285/H407.", "s": "Oxford; Sci Rep 2018"},
 "five_orient":   {"t": "A single PXR ligand can adopt near-isoenergetic modes — SR12813 was seen in five different orientations in the cavity.", "s": "Sci Rep 2018"},
 "five_subsites": {"t": "The pocket has five distinct subsites (Right/Left/Up/Down/Center) ligands extend into by size/shape — near-isoenergetic multi-subsite binding.", "s": "PMC"},
 "pocket_flex":   {"t": "PXR uses ligand-binding-domain flexibility as its core recognition mechanism; hyperforin expands the pocket by 250 A^3 — the pocket reshapes per ligand.", "s": "Biochemistry 2003"},
 "his407_hbond":  {"t": "His407 is a directional anchor that accepts a hydrogen bond from agonist hydroxyls (T0901317, SR12813).", "s": "Structure 2023"},
 "cs_weak":       {"t": "Thiocarbonyl sulfur (C=S) is a weak, anisotropic H-bond acceptor — ~20% weaker than a carbonyl O and accepting only along specific lone-pair directions.", "s": "ACS; PMC"},
 "af3_orient":    {"t": "AF3's dominant docking failure is ligand-pocket ORIENTATION error (right pocket, wrong pose): 58.6% of its errors on unseen ligands.", "s": "AF3 docking-failure analysis 2025"},
 "af3_hydrophobic":{"t": "AF3 is biased to place polar head-groups at polar residues and actively MIS-places hydrophobic moieties, misaligning the hydrophobic tail.", "s": "AF3 docking-failure analysis 2025"},
 "chai_degrade":  {"t": "Co-folding degrades sharply as the pocket diverges from training — Chai-1's median pose error rises ~2.4→5.7 A.", "s": "arXiv"},
}

# Campaign findings (our own results; traceable to memory memos in the repo's MEMORY index)
FINDINGS = {
 "selection_wall":   {"t": "Across the pool, best-of-N selection + cross-model z-hybrid is the ONLY thing that helped (0.46→0.50); no confidence signal or GT-trained selector reliably finds the good pose already in the pool.", "s": "project_what_moves_lddtpli; project_geometric_signals_refuted"},
 "plddt_noise":      {"t": "pLDDT/confidence is near-useless for pose selection here — it correlates ~0 with per-pose RMSD; the ~0.6 pLDDT-weight ceiling is a selection artefact, not a generation limit.", "s": "project_omnibus_iter4_final; project_signal_correlations_v1"},
 "size_bottleneck":  {"t": "Pose accuracy anti-correlates with size/flexibility (rho≈−0.55): fragments score ~0.73, large flexible analogs ~0.46.", "s": "project_size_is_the_bottleneck"},
 "input_fidelity":   {"t": "Input fidelity (native per-target sequence + deep MSA + sampling depth) dominates model choice: generic-sequence co-fold plateaus ~8 A while a native Boltz setup reaches 0.59 A on the same target.", "s": "project_posthoc_cross_model_gt"},
 "protenix_rescue":  {"t": "Swapping Protenix-v2 onto only the deepest failure-tail ligands (8 of 184) was the single best lever we found (best submission prot_rescue8=0.564).", "s": "project_combinatorial_ladder"},
 "templating_lost":  {"t": "Crystal-templating, pocket-anchor priors, and naive ensembles all LOST vs best-of-N on the GT harness — anchor priors mislead on these weak binders.", "s": "project_geometric_signals_refuted; project_activity_structure_coupling"},
 "activity_weak":    {"t": "82/111 drug-like analogs are the activity-track set and all 24 measured are PXR-INACTIVE — the whole 184 is a weak-binder co-folding problem, so anchor-based priors mislead.", "s": "project_activity_structure_coupling"},
}

# retrospective levers (what could have reached the pose) — id -> text, keyed by tags
RETRO = {
 "more_sampling":   "Generate far more poses (best-of-N raised the pool oracle) — but this only helps if paired with a selector, since the good pose is usually already in the pool.",
 "better_selector": "A better SCORING/selection function is the highest-value fix here: the pool typically already contains a near-native pose that confidence cannot pick.",
 "anchor_restraint":"A pocket-anchor restraint (Ser247/Gln285/His407) COULD bias the one polar handle — but we found anchor priors net-LOST on these weak binders, so use only as a soft, ligand-specific term.",
 "no_anchor_help":  "Anchor/template restraints are unlikely to help — there is no strong directional contact to restrain; the miss is hydrophobic/degenerate.",
 "input_fidelity":  "Higher input fidelity (native per-target MSA depth, more recycles) is the lever that moved poses most in our cross-model GT test.",
 "finetune":        "A PXR/NR-family fine-tune or a targeted model swap (Protenix-v2 on the failure tail) rescued the hardest cases better than a broader-but-shallower model sweep.",
 "cs_aware":        "A scoring term aware that C=S is a weak, directional acceptor could stop models over-rewarding a spurious S···H-bond.",
}

# ----------------------------------------------------------------------------- tag rules
# Percentile thresholds from the cohort
DIS_P75 = M["disagreement_mean_pw_rmsd"].quantile(.75)
DIS_P90 = M["disagreement_mean_pw_rmsd"].quantile(.90)
CONS_P25 = M["consensus_frac"].quantile(.25)


def polar_handles(r):
    return int(r["hbd"]) + int(r["hba"])


def tags_for(r, pose_ctx):
    """Return list of fired tags. Each tag: {id,label,why(data-grounded),claims,retro,sev}."""
    T = []
    mw, rotb, heavy, fsp3 = r["mw"], r["rotb"], r["heavy"], r["fsp3"]
    dis, cons = r["disagreement_mean_pw_rmsd"], r["consensus_frac"]
    ph = polar_handles(r)

    if rotb >= 5 or mw >= 400:
        T.append(dict(id="large_flexible", label="Large / flexible", sev=3,
            why=f"MW {mw:.0f}, {int(rotb)} rotatable bonds, {int(heavy)} heavy atoms. In a large uncharged pocket a flexible ligand has many near-isoenergetic placements, so the correct orientation is under-determined.",
            claims=["pocket_size","five_orient","five_subsites"], findings=["size_bottleneck"],
            retro=["more_sampling","better_selector"]))
    if ph <= 3 and fsp3 >= 0.3:
        T.append(dict(id="weak_anchor", label="Weak anchoring", sev=2,
            why=f"Only {int(r['hbd'])} H-bond donors + {int(r['hba'])} acceptors on a mostly aliphatic scaffold (Fsp3 {fsp3:.2f}). With so few polar handles for the Ser247/Gln285/His407 cluster, little pins the pose.",
            claims=["uncharged","anchor_cluster"], findings=["templating_lost"],
            retro=["no_anchor_help","better_selector"]))
    if fsp3 >= 0.45 or (r["clogp"] >= 3 and ph <= 3):
        T.append(dict(id="hydrophobic_misplace", label="Hydrophobic-tail risk", sev=2,
            why=f"A largely hydrophobic scaffold (Fsp3 {fsp3:.2f}, cLogP {r['clogp']:.1f}). AF3-class models bias polar groups toward polar residues and mis-place hydrophobic tails, exactly the group with no directional contact here.",
            claims=["af3_hydrophobic","uncharged"], findings=[],
            retro=["better_selector","input_fidelity"]))
    if int(r["nS"]) > 0:
        T.append(dict(id="sulfur", label="Sulfur H-bonding", sev=1,
            why=f"Contains {int(r['nS'])} sulfur atom(s). If any is a thiocarbonyl/thioether acceptor, models tend to mis-weight it — C=S is a weak, directional acceptor, unlike the carbonyl O they are trained on.",
            claims=["cs_weak"], findings=[],
            retro=["cs_aware"]))
    if dis >= DIS_P90:
        T.append(dict(id="very_high_disagreement", label="Models strongly disagree", sev=3,
            why=f"Cross-model pose disagreement {dis:.1f} A (top-decile of the 184). The 15 predictors place this ligand in genuinely different sub-sites — the pose is not pinned by the pocket.",
            claims=["five_subsites","pocket_flex","chai_degrade"], findings=["input_fidelity"],
            retro=["input_fidelity","finetune"]))
    elif dis >= DIS_P75:
        T.append(dict(id="high_disagreement", label="High model disagreement", sev=2,
            why=f"Cross-model pose disagreement {dis:.1f} A (upper quartile). Models diverge on where this ligand sits.",
            claims=["five_subsites","chai_degrade"], findings=["input_fidelity"],
            retro=["input_fidelity","better_selector"]))
    if cons <= CONS_P25:
        T.append(dict(id="low_consensus", label="No dominant pose", sev=2,
            why=f"Only {cons:.0%} of models fall in the largest pose cluster — there is no majority pose to trust as a consensus answer.",
            claims=["five_orient","pocket_flex"], findings=["selection_wall"],
            retro=["better_selector"]))
    if r["population"] == "fragment" and heavy <= 16:
        T.append(dict(id="fragment", label="Fragment — too few contacts", sev=2,
            why=f"A {int(heavy)}-heavy-atom PanDDA fragment. Too few contacts to pin a pose in a large promiscuous pocket, and the apo PanDDA site is dissimilar to the models' training complexes.",
            claims=["pocket_size","chai_degrade"], findings=["size_bottleneck"],
            retro=["more_sampling","better_selector"]))
    # selection decoupling (confidence disagrees with consensus)
    if pose_ctx.get("plddt_pick") and pose_ctx["plddt_pick"] != r["medoid_model"]:
        T.append(dict(id="selection_risk", label="Confidence ≠ consensus", sev=3,
            why=f"The pLDDT-preferred pose comes from {pose_ctx['plddt_pick']}, but the cross-model consensus (medoid) is {r['medoid_model']}. Confidence and consensus disagree — the hallmark of the selection wall (true in 172/184 ligands).",
            claims=["af3_orient"], findings=["plddt_noise","selection_wall"],
            retro=["better_selector"]))
    return T


def main():
    # per-structure pLDDT pick (which model's pose confidence would choose)
    picks = ML.loc[ML.groupby("structure")["lig_plddt"].idxmax()][["structure", "model"]] \
              .set_index("structure")["model"].to_dict()
    dis_rank = M["disagreement_mean_pw_rmsd"].rank(pct=True)

    out = {}
    for i, r in M.iterrows():
        sid = r["structure"]
        ctx = {"plddt_pick": picks.get(sid)}
        T = tags_for(r, ctx)
        # per-pose journey summary from long table
        sub = ML[ML.structure == sid].sort_values("rmsd_to_medoid")
        journey = [dict(model=x.model, rmsd_to_medoid=round(float(x.rmsd_to_medoid), 2),
                        lig_plddt=round(float(x.lig_plddt), 3), is_medoid=bool(x.is_medoid))
                   for x in sub.itertuples()]
        retro_ids, claim_ids, finding_ids = [], [], []
        for t in T:
            retro_ids += t.get("retro", []); claim_ids += t.get("claims", []); finding_ids += t.get("findings", [])
        retro_ids = list(dict.fromkeys(retro_ids)); claim_ids = list(dict.fromkeys(claim_ids)); finding_ids = list(dict.fromkeys(finding_ids))
        out[sid] = dict(
            sid=sid, name=r.get("mol_name") or sid, population=r["population"], series=r.get("series"),
            smiles=r["smiles"], murcko=r.get("murcko"),
            desc=dict(mw=round(float(r.mw),1), clogp=round(float(r.clogp),2), tpsa=round(float(r.tpsa),1),
                      hbd=int(r.hbd), hba=int(r.hba), rotb=int(r.rotb), heavy=int(r.heavy),
                      fsp3=round(float(r.fsp3),2), narom=int(r.narom), nrings=int(r.nrings), nS=int(r.nS)),
            activity=dict(pEC50=None if pd.isna(r.pEC50) else round(float(r.pEC50),2),
                          emax=None if pd.isna(r.emax) else round(float(r.emax),1)),
            xmodel=dict(disagreement=round(float(r.disagreement_mean_pw_rmsd),2),
                        disagreement_pctile=round(float(dis_rank[i]),2),
                        consensus_frac=round(float(r.consensus_frac),2),
                        medoid_model=r.medoid_model, plddt_pick=picks.get(sid),
                        selection_decoupled=bool(picks.get(sid) != r.medoid_model),
                        mean_lig_plddt=round(float(r.mean_lig_plddt),1), n_models=int(r.n_models)),
            tags=[{k:t[k] for k in ("id","label","sev","why","claims","findings")} for t in T],
            retrospective=[{"id":rid,"text":RETRO[rid]} for rid in retro_ids],
            journey=journey,
            has_crystal=False,
        )
    payload = dict(claims=CLAIMS, findings=FINDINGS,
                   thresholds=dict(dis_p75=round(float(DIS_P75),2), dis_p90=round(float(DIS_P90),2),
                                   cons_p25=round(float(CONS_P25),2)),
                   n=len(out), compounds=out)
    os.makedirs("data/processed/posthoc", exist_ok=True)
    json.dump(payload, open("data/processed/posthoc/compound_analysis.json","w"))
    # quick provenance summary
    from collections import Counter
    tagc = Counter(t["id"] for c in out.values() for t in c["tags"])
    print(f"wrote compound_analysis.json: {len(out)} compounds")
    print("selection-decoupled:", sum(c['xmodel']['selection_decoupled'] for c in out.values()), "/", len(out))
    print("tag frequencies:", dict(tagc.most_common()))
    print("compounds with 0 tags:", sum(1 for c in out.values() if not c['tags']))


if __name__ == "__main__":
    main()
