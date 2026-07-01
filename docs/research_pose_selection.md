# PXR Pose-Selection Research Synthesis (2026-06-23)

Deep-web-search synthesis (research agent, 10 tool-uses, ~4-5 cycles). Question: best custom scoring/selection
function to pick the highest-LDDT-PLI pose from a large cofolded pool for the PXR challenge.

## Bottom line
2024-2026 literature is unanimous and *bad news for a clever global reranker*: native confidence ranking and
cross-model **consensus largely do NOT beat random** on cofolding pose pools, and consensus can be **actively
harmful** (agreeing models share correlated errors). Our harness result (consensus≈random, oracle≫realized) is
the EXPECTED result. Proven headroom = (a) **interface-restricted confidence** fields, (b) **more sampling**
to raise the oracle, (c) **targeted per-ligand rescue**. No public ready-to-use self-supervised LDDT-PLI
estimator exists to drop in.

## Ranked recommendations (for our setup)
1. **Interface-RESTRICTED confidence re-rank — HIGHEST EV, ~1 day, no GPU.** We likely select on GLOBAL
   ranking_score/iptm/complex_plddt which DILUTE the small ligand interface across the large LBD. Switch to:
   - AF3: **ligand-local pLDDT** (best single signal, AUC~0.76) + **chain-pair-restricted PAE / ipSAE** over
     the protein↔ligand block. (ranking_score=0.8 ipTM+0.2 pTM+... is the WRONG field for pose quality.)
   - Boltz-2: ligand_iptm + complex_ipde (we ALREADY use IPDE — correct, keep). NOT Boltz-2 affinity (r=-0.03).
   - Impl: parse PAE JSONs; compute chain-pair PAE MANUALLY over protein↔ligand block. Stock DunbrackLab/IPSAE
     ipsae.py EXCLUDES ligands (label_seq_id=".") -> must ADAPT it.
2. **Sampling->oracle->failure-tail rescue — proven, already in motion.** top-confidence pose ≈ random;
   best@1 79%->best@5 86%->best@20 89% (Runs N' Poses). Our 0.5629/0.5640 targeted Protenix rescue is the
   literature-correct pattern. Use cross-model diversity to WIDEN the pool, NOT to vote.
3. **PoseBusters as a validity GATE (not ranker) — cheap insurance, ~half day.** `pip install posebusters`;
   ~20 RDKit checks (valence/stereo/clash). Drop physically-broken poses BEFORE ranking + QC final 184 PDBs
   (ties into our scoring-fail work). Relaxation can't fix chirality -> gate, don't fix.
4. **RMSD-Pred (eightmm/RMSD-Pred, JCIM 2025) — model-agnostic rescorer, GT-gated, ~1-2 days.** GatedGCN:
   protein PDB + multi-pose SDF -> predicted RMSD + P(RMSD>2A) per pose. Bundled weights; HAS OOD validation
   (CleanSplit gap <2%). Use RMSD-Pred for POSES; do NOT use its BA-Pred affinity. Validate on 35-53 holo GT
   with PROTEIN-CLUSTERED CV first (CASF-trained backbone -> off-domain risk on flexible analogs).
5. **gnina CNNscore — orthogonal tiebreaker, leave-protein-out gated, ~1 day.** `gnina --score_only
   --cnn_scoring rescore`; CNNscore=P(pose<=2A). Best combined with RMSD-Pred via RANK-SUM (each weak alone;
   win is in orthogonal-scorer consensus). Use CNNscore not CNNaffinity.

## PXR-specific term (SOFT tiebreaker, analogs ONLY — additive bonus, never penalize absence)
- Primary: ligand H-bond donor/acceptor within ~3.5A of **Ser247 Oγ** (α3, agonist-defining) and/or **Gln285**.
- Secondary: **His407** (α10, flexible). Aromatic burial: ring centroid within ~4.5-5.5A of Phe288/Trp299/Tyr306.
- CAUTION: LBP huge (>1600 A³), promiscuous; **PanDDA FRAGMENTS often engage ZERO canonical anchors** -> this
  term MISLEADS on fragments (matches our prior anchor-prior finding). Helps mostly drug-like analogs.

## One careful GT-gated shot
ProLIF directional interaction-fingerprint recovery vs closest holo-analog reference — the only signal shown
genuinely ORTHOGONAL to RMSD/confidence (Errington 2024). Directional bits only, Tanimoto to reference. Helps
fragments WITH a holo reference; as a reference-free selector it's unproven.

## AVOID (evidence-backed)
Global ipTM/pLDDT/ranking_score for pose ranking; Boltz-2 affinity; **consensus voting/clustering/centroid**
(correlated-error trap); ligand-only strain/MMFF min (bound conformers legitimately strained — do protein-
restrained COMPLEX relax if any); many-feature learned selectors trained on the 35-53 GT (overfit).

## Overfitting reality check
Selector AUCs 0.73-0.76 need HUNDREDS of systems; on 35-53 holos any ±0.05 "win" is NOISE (our pocket_plddt
rose only +0.0015 on 53-holo — the cautionary case). Use scaffold/protein-clustered CV; report gap-to-oracle
and gap-to-random; **trust the leaderboard over the tiny GT.**

## Suggested concrete pipeline (next few days)
Per ligand: (1) PoseBusters-gate pool -> (2) re-rank survivors by **ligand-local pLDDT / chain-pair-PAE
(ipSAE-adapted) + existing Boltz IPDE** -> (3) add **RMSD-Pred ∧ gnina CNNscore** rank-sum (orthogonal) ->
(4) optional small PXR-anchor bonus on ANALOGS only. Validate the ORDERING on 35-53 holo GT (protein-clustered
CV) before spending a ladder slot.

## Key sources
ipSAE PMC11844409 · DunbrackLab/IPSAE · Pearl arXiv 2510.24670 · Runs N' Poses (Nat Struct Mol Biol 2026) ·
Mac1 cofolding eval eLife 110475 / PMC12776374 · RMSD-Pred github.com/eightmm/RMSD-Pred / PMC13080981 ·
gnina github.com/gnina/gnina / v1.3 JCheminf · PoseBusters RSC d3sc04185a / PyPI · ProLIF arXiv 2409.20227 ·
GenScore github.com/sc8668/GenScore · PXR pocket PMC9563780 / PMC8864553 / Science promiscuity (1060762) ·
Inductive Bio gnina rescore blog.
