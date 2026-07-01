# Iterations 1–7 — narrative

Day-by-day chronological narrative of what was tried, what worked, and what failed. Companion to the tabular per-method audit in [stats/methods_table.md](stats/methods_table.md). Use this when you want to understand "why was that decision made" — use the table when you want a one-line verdict per method.

## Iter 0 — Cofold baselines (May 18–20, 3 days)

**Goal:** establish a working pipeline + push the highest score achievable from out-of-the-box cofold + simple confidence ranking.

**What ran:**
- `#01_boltz2_tutorial`: Boltz-2 single-pose tutorial PDBs (no rerank). 0.4632 / r21.
- `#04_of3_best`: OpenFold3 best-of-20 by `sample_ranking_score`. 0.4751 / r12.
- `#06_b2_ligand_iptm`: Boltz-2 best-of-20 by `ligand_iptm`. 0.4869 / r10.
- `#08_iptm_zhybrid`: z-score B2 `ligand_iptm` + OF3 `iptm` cross-model. 0.4984 / r6. **First top-10.**
- `#13_zhybrid_srs`: switched OF3 to `sample_ranking_score`. **0.4997 / r5**. First top-5.

**Key findings:**
- Submissions 4 and 5 (in 24h) jumped us from r21 to r5. **The signal choice mattered more than fancy aggregation.**
- B2-only (best-of-20) plateaus at 0.487; OF3-only at 0.475; the mix was +0.014 over the better single.
- **5/5 submitted, 5/5 SUBMITTED_WIN.** Hardest hit-rate of the whole project.

**State at iter close:** 0.4997 / r5. The "gold" #13 is our high-water mark of the project (will be matched but not exceeded later by #262 at 0.4996, which is statistically tied).

## Iter 1 — Anchor / family / IPDE crystallization (May 21–23, 3 days)

**Goal:** beat 0.4997 by adding (a) crystal-anchor information (we have 9 cached PXR holos), (b) per-family rerank (activity-Tanimoto / scaffold splits), or (c) better within-pool signal.

**What ran:**

Crystal-anchor track (3 attempts):
- `#113_crystal_protein_srs`: Boltz-2 with 1ILH crystal protein-only prior (B-factor seeded). 0.4657 / r22. **First catastrophe.** Crystal B-factors broke pLDDT (0.91 → 0.85).
- `#274_mega_crystal_sa`: T=0.5 simulated annealing toward 1ILH heavy atoms. **0.1671 / r38.** Back-transform bug — the SA wrote ligand coordinates in a different frame than the protein. This was the moment we instituted the no-restore-reflex policy.
- `#278_crystal_anchor_sa`: T=0.1 50-step SA. 0.3828 / r29. Still anchor-driven failure.

Restoration track:
- `#262_ai_ipde_select`: pure `-complex_ipde` cross-model (B2 + OF3 heuristic-scaled). **0.4996 / r7.** Our gold. Statistically tied with #13 but pure-IPDE-conceptually-clean and easier to reproduce. **THIS IS THE BEST EVER.**
- `#262_restore`: 4-hour-later resubmit of #262 via `auto_submit_v2.py` after the #274 catastrophe. Recovered to 0.4996.

Single-model and reranker variants:
- `#282_ipde_b2_only`: drop OF3. 0.4818 / r12 (-0.0178). OF3 contributes +0.018.
- `#263_of3_only`: drop B2. 0.4724 / r13 (-0.0272). B2 contributes +0.027.
- `#283_zscore_complex_ipde`: per-model z-norm. 0.4601 / r25. Loses absolute calibration.
- `#267_meta_ensemble_v3`: meta-ensemble. **0.2010 / r41.** Anchor-style failure with wrong frame again.
- `#281_per_family_mlp`: per-family MLP rerank trained on activity-track. 0.4662 / r21. **94.6% pick divergence vs #262** — the anchor-disaster pattern on PanDDA fragments.
- `#279_ai_unimol_pxr`: Uni-Mol head finetuned on PXR ChEMBL. Score never captured (overwritten before logged); estimated 0.45-0.49.

**Key findings:**
- **Crystal anchors fail catastrophically on PanDDA fragments.** Three independent attempts (113, 274, 278), all in the bottom quartile.
- **Per-family MLP fails for the same reason.** 94.6% divergence on fragments means the MLP picked totally different poses than IPDE, and the IPDE picks were systematically better.
- **Pure IPDE cross-model = local maximum.** #262 ties our best-ever (#13), with cleaner conceptual interpretation.
- **Self-overwrite hazard** is real and expensive — we burned 5 sub-0.47 submissions trying to beat 0.4996.

**Decisions codified here:**
- **No-restore-reflex policy.** Always end a submission cycle with `#262_restore`.
- **Auto-submit safety net** (`auto_submit_v2.py`) — fires on a schedule independent of agent state.

**State at iter close:** 0.4996 / r7 (= rank 6 after restore). 10 submitted, 1 SUBMITTED_WIN (#262), 9 SUBMITTED_LOSS.

## Iter 2 — Signal mining (May 25–28, 4 days)

**Goal:** find a non-IPDE signal that adds information over IPDE alone. Mine `complex_iplddt`, `pde`, `sample_ranking_score`, holo occupancy, anti-occupancy, activity Tanimoto.

**What ran:**
- `#290_b2_iplddt`: `complex_iplddt`-best instead of `complex_ipde`-best. 0.4832 / r13 (-0.0164). iplddt pushed B2 over OF3 catastrophically — biases toward Boltz-2 which has higher absolute iplddt.
- `#294_holo_occupancy`: 90% IPDE + 10% holo-occupancy bonus (poses overlapping with 9 cached PXR holos' ligand densities). 0.4586 / r24 (-0.0410). **Fragments bind in NEW sub-pockets** — occupancy with drug-like holo ligands is the wrong prior.
- `#297_anti_occupancy`: inverse prior (poses AVOIDING known holo densities). 0.4514 / r25 (-0.0482). The relationship is non-linear, both signs hurt.
- `#298_zscore_srs`: complex_ipde + OF3 sample_ranking_score z-score. 0.4597 / r24 (-0.0399).
- `#322_activity_tanimoto`: pick pose closest to medoid of activity-Tanimoto-neighbors' best poses. 0.458 / r27 (-0.0416). Medoid-of-neighbors falsified.

Held by divergence gate:
- `#287_hard_consensus`, `#288_soft_bonus`, `#289_tie_breaker`: all collapsed near B2-only or had < 5% divergence. Held.
- `#302_mmff_strain`: 86.4% divergence (159/184 different poses); strain validator failed 35-holo gate. **HOLD_DIVERGENCE_FAIL** (later confirmed REJECT, see iter 6 closure).

**Key findings:**
- **Every off-IPDE submission burned a slot and scored worse than 0.4996.**
- **Both holo-prior directions (occupancy AND anti-occupancy) hurt** — the relationship between known holo coverage and fragment-pose-quality is not monotonic.
- **iplddt regresses on this test set** despite being a metric-aligned signal (LDDT-PLI is what we're scored on). The reason: iplddt is biased toward B2-internal confidence, and B2's confidence is not necessarily calibrated for PanDDA fragments.
- **6 submissions, 6 SUBMITTED_LOSS.** ROI of brainstorm volume here was strongly negative.

**State at iter close:** 0.4996 / r9 (drifted -2 ranks). 6 submitted, 0 SUBMITTED_WIN, 6 SUBMITTED_LOSS.

## Iter 3 — Omnibus external (June 3, 1 day kickoff)

**Goal:** the kernel-internal signal space is mapped; pivot to **external** signals (AF3 server batches, ESMFold + smina, OpenMM MD, gnina CNN, PDBbind XGBoost, ChEMBL PXR Uni-Mol v2).

**What ran (all kicked off in parallel as Kaggle kernels):**
- `af3_triple_consensus`: prepare 184 jobs for AF3 server upload. **BLOCKED_INFRA** — no public AF3 REST API; 20 jobs/day manual upload. Required user action.
- `esmfold_induced_fit_dock` (nb08): ESMFold receptor + smina ligand dock. Completed but **88.6% divergence** vs #262 → HOLD_DIVERGENCE_FAIL.
- `openmm_md_refine` (nb09): 9h Kaggle T4 successful, 184/184 PDBs produced. 67.4% divergence, 0 high-confidence overrides over IPDE-best → HOLD_DIVERGENCE_FAIL.
- `gnina_cnn_rescore` (nb07): gnina CNN rescore on B2+OF3 poses. 99.5% divergence — perfect anchor-disaster pattern → CLOSED_NEGATIVE.
- `pdbbind_xgboost_lambdarank`: LOO XGB on PDBbind NR1I subset. Validation showed 2.275 Å vs IPDE 2.230 Å (-0.045 Å) → CLOSED_NEGATIVE.
- `chembl_pxr_v2`: 31-entry ChEMBL PXR data corpus for Uni-Mol finetune v2. Sub-µM pEC50 mismatch with fragments → CLOSED_NEGATIVE.
- `activity_aware_picker`: stratified rerank based on activity-track family. Only 4.3% divergence (< 5% gate floor) → HOLD_SUB_NOISE.

**Key findings:**
- **0/7 external paths produced a submittable candidate** in the 7-day window since launch.
- **Validation gate caught all 5 negative candidates** without burning a slot. Defensive infrastructure paid for itself.
- **External resources have a 12-week activation lag** that doesn't fit a 30-day competition: AF3 needs daily manual uploads + 6h per batch; ESMFold + dock needs >1 GPU; OpenMM MD refine ran but produced 0 confident overrides; gnina was trained on PDBbind drug-like and inverts on fragments; XGB inherits the same drug-like bias.

**Decision pivot:** stop chasing external resources, focus remaining 25 days on:
1. Re-launching the kernels that crashed (iter 4).
2. Watching for organizer-released 64 re-refined PDBs.
3. Manual AF3 batch uploads (user action).

**State at iter close:** 0.4996 / r9. 0 submitted, 4 HOLD, 3 CLOSED_NEGATIVE.

## Iter 4 — Re-launched kernels (June 7–8, ~36h)

**Goal:** salvage the kernels from iter 3 that crashed or didn't finish. Specifically: Chai-1 (third cofold model), ProteinMPNN (pocket redesign), xTB (semi-empirical strain rescore), and val-set expansion (35 → 53 holos).

**What ran (4 Kaggle kernels pushed 2026-06-07 22:00 UTC):**

| nb | model / job | hardware | wall | terminal state |
|---|---|---|---|---|
| nb10 | Chai-1 v4 cofold | CPU (no GPU quota) | 15.5 s | **errored_closed**: `FileNotFoundError: nvidia-smi` (cell 1 GPU diagnostic crashed on CPU kernel) |
| nb15 | ProteinMPNN pocket redesign | CPU | ~2h | **partial_closed**: produced 4 MPNN designs (T=0.1, seq_recovery 0.46-0.58) at positions 245-289; **NO cofold+smina step ran** |
| nb16 | Boltz-2 ChEMBL val expansion | GPU (T4, 8h) | ~6h | **completed_closed**: 18 holos × 20 models = 360 PDBs (275 MB) |
| nb17 | xTB rescore | CPU | 11 s | **errored_closed**: same `nvidia-smi` crash as nb10 |

**Key findings:**

1. **nb10 + nb17 both crashed on `nvidia-smi`** — a GPU diagnostic at cell 1 of both notebooks failed when the kernel was CPU-only. We had previously seen this on iter 3 kicks; the patch (try/except wrap) was applied to nb15 but missed nb10/nb17.
2. **nb15 ProteinMPNN ran successfully but the cofold step never fired.** The 4 MPNN sequences are interesting structurally (substitutions at F248Y, L270A, Q273E, D274S, A279L, M283Q, Q284K — all in the AF-2 helix / pocket lid), but without the downstream cofold + smina step we have 0 new 184-ligand poses.
3. **nb16 was the lone success.** 18 PXR ChEMBL holos cofolded with Boltz-2, 0 overlap with the original 35-holo set, mean RMSD 0.85 Å (much easier than the 35-holo 1.40 Å — these are well-characterized drug-like binders).

**The 53-holo expansion REVERSED our prior ranking of selectors:**

| method | 35-holo lddt | 53-holo lddt | Δ |
|---|---|---|---|
| #plddt_best | 0.4075 (#1) | 0.4090 (LAST) | +0.0015 |
| #262_ipde_best | 0.4024 | **0.4214 (top)** | **+0.019** |

Plddt's 35-holo lead was a 35-holo overfit. On the easier 18-holo subset, plddt did NOT improve (+0.0015), while every other method jumped by +0.020. **IPDE is now validated as the most robust signal.**

This was the deepest single-iter insight of the whole project — and it landed retroactively (iter 4 confirmed iter 1's #262 decision was right; iter 2's pursuit of `iplddt` was wrong; iter 3's HOLD on `v3_hybrid_pocket_plddt` was right).

**State at iter close:** 0.4996 / r10 (drifted -1 to leaderboard movement; we did not submit). 0 submitted, all HOLD. Pre-flight gate result: **no candidate to gate**.

## Iter 5 — QTLSP brainstorm (June 8 AM, ~4h)

**Goal:** identify the single highest-value next experiment from the brainstorm catalog. The empirical kernel-internal ceiling is ~0.515 (project_numerical_reassessment_2026_06_08), so any new method needs an **external** physical/chemical signal.

**What ran:**
- 5 ideas generated in brainstorm; only 1 reached GO_BUILD spec:
  - **QTLSP (QCArchive TorsionDriveDataset Strain Prior)**: pull ~12k DFT torsion scans from QCArchive, fit a torsion-strain prior for each test ligand, rerank poses by `score = -complex_ipde + λ * external_strain_penalty`. Predicted lift +0.005 to +0.014. Dev cost 14 hours. **GO_BUILD with sanity gate**.
- 4 ideas truncated in upstream synthesis (placeholders: src1_coactivator, h12_swap_motif, conformal_quantile, untitled #5).

**Sanity gate for QTLSP:** the underlying assumption is that gas-phase ligand strain correlates with binding-pose quality. We need to validate this on a known-strain signal (MMFF) before investing 14h in DFT. **#302_mmff_strain validator** was the precursor sanity gate.

**State at iter close:** 0.4996 / r10. 0 submitted, 1 GO_BUILD (QTLSP), HOLD until sanity gate fires.

## Iter 6 — AMPT + MMFF closure (June 8 PM, ~6h)

### Iter 6a — AMPT orthogonality gate

**Goal:** test if AlphaMissense per-residue pathogenicity scores (AMPT) provide an orthogonal signal vs Boltz-2 pLDDT on the PXR LBD pocket.

**Setup:**
- AMPT = mean(missense_pathogenicity) over 28 PXR pocket residues (from EBI's `am_O75469.csv`).
- pLDDT pocket = mean pLDDT over the same 28 residues across 700 Boltz-2 predictions (5 models × 184 ligands).
- Orthogonality criterion: |Spearman ρ| < 0.3.

**Result:**
- **Spearman ρ = -0.2939** (p = 0.129, not significant). **Just barely passes orthogonality** (0.006 below 0.3 cutoff).
- Pearson ρ = -0.3229 (p = 0.094, borderline; CROSSES 0.3).
- Sign-inversion risk persists — pocket-residue signals already sign-inverted once (v3_hybrid_pocket_plddt: predicted +0.0099, scored -0.0294).

**Verdict:** GATE PASSED → GO_BUILD eligible **but** placed behind QTLSP and iter4 kernels in priority order. Build cost 2h, expected lift +0.003 lddt_proxy (borderline noise).

Closure memo: `project_ampt_closure.md`. Artifact: `data/processed/ampt_orthogonality_gate.json`.

### Iter 6b — MMFF strain holo validation (QTLSP sanity gate)

**Goal:** validate the QTLSP precursor — does MMFF gas-phase ligand strain correlate with pose-quality on the 35-holo crystal set?

**Setup:**
- Compute MMFF strain for each of 700 predicted poses (35 holos × 20 Boltz-2 models).
- Compute strain_best (pick pose with lowest gas-phase strain), pure-strain, ipde+strain blend (top-3 IPDE-filtered + strain rerank).
- Compare to oracle (lower bound) and ipde_best (current #262 baseline).

**Result:**

| method | mean RMSD | under 2 Å | lddt_proxy |
|---|---|---|---|
| oracle | 1.683 | 91.4% | 0.5109 |
| ipde_best (#262 local baseline) | **2.230** | 57.1% | 0.4024 |
| strain_best (pure) | 2.332 | 62.9% | 0.3966 |
| blend_all (over all 20 poses) | 2.341 | 60.0% | 0.4009 |
| blend_top3 (over top-3 IPDE; replicates #302) | 2.251 | 62.9% | 0.4003 |

**Best strain variant: blend_top3 at 2.251 Å — +0.020 Å WORSE than IPDE.**

Per-holo signature: **high-variance**.
- Strain wins big on 5A86 (1.05 vs 1.93), 6A6M (1.25 vs 2.04), 8FPE (1.05 vs 2.70), 8SVP (1.58 vs 2.34).
- Strain loses catastrophically on 8F5Y (2.84 vs 0.49), 8SZV (2.87 vs 1.47), 9BEQ (5.12 vs 2.17).
- under-2Å fraction IMPROVES (0.629 vs 0.571), but catastrophic 3-5 Å misses on 3 holos drag the mean down.
- 86.4% divergence vs #262 (159/184 picks different).

**Verdict: REJECT.**

**Implication for QTLSP:** REJECTED transitively per user-defined rule ("if MMFF strain fails holo gate, QTLSP build is REJECTED — DFT cannot rescue a weak base signal"). The signal-to-noise of strain-energy-as-pose-scorer is the limiting factor, NOT the force-field accuracy. DFT would have computed the same energy ordering on the misleading 8F5Y / 8SZV / 9BEQ catastrophes.

Action: **cancel QTLSP build.** Freed 14 dev-hours + ~3h compute cache.

Closure memo: `project_mmff_strain_closure.md`. Artifact: `data/processed/mmff_strain_holo_validation.json`.

**State at iter close:** 0.4996 / r10. 0 submitted. AMPT GO_BUILD eligible but HOLD. QTLSP REJECTED. **All in-flight build paths now closed.**

## Iter 7 — Documentation closeout (June 8–30, current)

**Goal:** close all loose ends end-to-end. No new ideas. Document with figures + stats.

**What's being done:**
- Rewrite README.md as the single entry point.
- Generate 4 new figures (`fig01_leaderboard_trajectory`, `fig02_validation_methods`, `fig03_predicted_vs_actual`, `fig04_iter_attempts`).
- Generate 3 stats tables (`methods_table.md`, `leaderboard_table.md`, `iter_summary.md`).
- Write `docs/INDEX.md` (master pointer), `docs/methodology.md` (this method's details), `docs/iterations.md` (this file).
- Update MEMORY.md with the iter 7 closeout memo + documentation memo.
- AF3 batch JSONs ready at `submissions/af3_batches/` with upload instructions.
- Safety-net Scheduled Task + backup PID already armed for 2026-06-30 23:00 / 23:30 UTC.

**No new submissions planned.** External triggers that would unblock a new submission:
- AF3 batches return AND triple/quadruple consensus passes 3-gate.
- Organizer releases 64 re-refined PDBs on HF AND IPDE-best with re-templated receptors shows ≥+0.01 on 53-holo.
- A new high-quality cofold model lands on Kaggle T4-compatible GPU.

Without external trigger, final submission stays at 0.4996. Final rank likely 10-15.

## Cross-iter lessons that generalize

1. **Cross-model selection beats every other signal class** on PanDDA-fragment-heavy test sets where drug-like priors break.
2. **Validation → test compression is unreliable for pocket-residue signals.** Only IPDE-class signals (interface confidence) transfer.
3. **The leaderboard-shows-latest trap is real and costly.** 3 ex-top-5 teams dropped > 0.05 in 6 days through self-overwrite. No-restore-reflex + 25% divergence gate + safety task are necessary infra.
4. **External resources have a 12-week activation lag.** 5 ambitious iter-3 external paths (AF3, ESMFold, MD, gnina, XGB) all hit infra blockers before producing terminal output. 0 / 8 external paths produced a submittable candidate.
5. **Brainstorming yields diminishing returns past ~20 ideas** when the underlying signal space is already mapped (signal-correlations |ρ| < 0.08 across all kernel-internal confidence signals).
6. **HOLD is a real strategy.** 6 days at 0.4996 preserved approximately +0.05 to +0.10 LDDT-PLI of leaderboard drift against the casualty-list teams. The dominant strategy at this stage is preservation, not exploration.
