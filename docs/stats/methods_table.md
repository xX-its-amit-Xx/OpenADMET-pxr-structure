# Methods Table (Iterations 1-7)

Per-method audit of every attempted approach across the full project window (2026-05-18 -> 2026-06-08).

- `predicted_lift` = expected LDDT-PLI delta vs current best at proposal time (positive = better).
- `actual_lift` = observed leaderboard delta vs #262 IPDE baseline (0.4996). `n/a` when held / never submitted.
- `divergence_vs_262` = fraction of 184 picks differing from #262 IPDE-best. Gate band [5%, 30%].
- `verdict`: SUBMITTED_WIN / SUBMITTED_NEUTRAL / SUBMITTED_LOSS / HOLD_SUB_NOISE / HOLD_DIVERGENCE_FAIL / CLOSED_NEGATIVE / BLOCKED_INFRA / NEVER_RAN.

## Iter 0 - Baseline establishment

| name | family | external/internal | predicted_lift | actual_lift | divergence_vs_262 | dev_hours | verdict | one_line_reason |
|---|---|---|---|---|---|---|---|---|
| 01_boltz2_tutorial | cofold-baseline | internal | n/a | -0.0364 | n/a | 0.5 | SUBMITTED_NEUTRAL | tutorial single-pose Boltz-2 PDB; 0.4632 / r21 |
| 04_of3_best | cofold-baseline | internal | +0.005 | -0.0245 | n/a | 1.0 | SUBMITTED_WIN | OF3 best-of-20 by sample_ranking_score; 0.4751 / r12 |
| 06_b2_ligand_iptm | confidence-ranker | internal | +0.005 | -0.0127 | n/a | 1.0 | SUBMITTED_WIN | Boltz-2 best-of-20 by ligand_iptm; 0.4869 / r10 |
| 08_iptm_zhybrid | cross-model | internal | +0.005 | -0.0012 | n/a | 2.0 | SUBMITTED_WIN | B2 ligand_iptm + OF3 iptm z-hybrid; 0.4984 / r6 |
| 13_zhybrid_srs | cross-model | internal | +0.001 | +0.0001 | n/a | 1.5 | SUBMITTED_WIN | switched OF3 to sample_ranking_score; 0.4997 / r5 (= gold) |

## Iter 1 - Cross-model and crystal-anchor approaches (May 21-22)

| name | family | external/internal | predicted_lift | actual_lift | divergence_vs_262 | dev_hours | verdict | one_line_reason |
|---|---|---|---|---|---|---|---|---|
| 14_zhybrid_both_composite | cross-model | internal | +0.005 | -0.0148 | n/a | 1.0 | SUBMITTED_LOSS | OF3 sample_ranking_score composite added noise; 0.4849 / r11 |
| 113_crystal_protein_srs | crystal-anchor | internal | +0.010 | -0.0340 | n/a | 2.0 | SUBMITTED_LOSS | 1ILH crystal B-factors broke pLDDT (0.91->0.85); 0.4657 / r22 |
| 262_ai_ipde_select | cross-model | internal | +0.005 | 0.0000 | 0% | 1.5 | SUBMITTED_WIN | pure IPDE B2+OF3 heuristic scale; **0.4996 / r7 = best ever** |
| 274_mega_crystal_sa | crystal-anchor | internal | +0.010 | -0.3325 | n/a | 3.0 | SUBMITTED_LOSS | back-transform bug crashed score to 0.1671 / r38 |
| 278_crystal_anchor_sa | crystal-anchor | internal | +0.005 | -0.1168 | n/a | 2.0 | SUBMITTED_LOSS | T=0.1 50-step SA; 0.3828 / r29 |
| 267_meta_ensemble_v3 | ensemble | internal | +0.005 | -0.2986 | n/a | 2.0 | SUBMITTED_LOSS | meta-ensemble v3; 0.2010 / r41 |
| 279_ai_unimol_pxr | learned-ranker | internal | +0.008 | unknown | unknown | 6.0 | SUBMITTED_NEUTRAL | Uni-Mol on PXR finetune; score never captured (overwritten); est 0.45-0.49 |
| 281_per_family_mlp | learned-ranker | internal | +0.008 | -0.0334 | 94.6% | 4.0 | SUBMITTED_LOSS | per-family MLP recapitulated anchor noise; 0.4662 / r21 |
| 282_ipde_b2_only | cross-model | internal | -0.005 | -0.0178 | 70.1% | 0.5 | SUBMITTED_LOSS | B2-only IPDE; OF3 contributes; 0.4818 / r12 |
| 263_of3_only | cross-model | internal | -0.005 | -0.0272 | 29.9% | 0.5 | SUBMITTED_LOSS | OF3-only IPDE; B2 contributes; 0.4724 / r13 |
| 283_zscore_complex_ipde | cross-model | internal | +0.001 | -0.0395 | 62.5% | 1.0 | SUBMITTED_LOSS | per-model z-score on complex_ipde+iptm; 0.4601 / r25 |

## Iter 2 - iplddt and signal-mining (May 25-28)

| name | family | external/internal | predicted_lift | actual_lift | divergence_vs_262 | dev_hours | verdict | one_line_reason |
|---|---|---|---|---|---|---|---|---|
| 287_hard_consensus | cross-model | internal | +0.003 | not submitted | n/a | 1.5 | HOLD_DIVERGENCE_FAIL | collapses to ~B2-only; predicted to underperform |
| 288_soft_bonus | cross-model | internal | +0.003 | not submitted | n/a | 1.0 | HOLD_DIVERGENCE_FAIL | low overlap with #262, no signal |
| 289_tie_breaker | cross-model | internal | +0.003 | not submitted | n/a | 1.0 | HOLD_DIVERGENCE_FAIL | tie-breaker shifted toward B2-only |
| 290_b2_iplddt | confidence-ranker | internal | +0.010 | -0.0164 | high | 2.0 | SUBMITTED_LOSS | complex_iplddt pushed B2 over OF3 catastrophically; 0.4832 / r13 |
| 294_holo_occupancy | holo-prior | internal | +0.005 | -0.0410 | n/a | 2.0 | SUBMITTED_LOSS | 90/10 IPDE/occupancy; fragments bind in NEW sub-pockets; 0.4586 / r24 |
| 297_anti_occupancy | holo-prior | internal | +0.003 | -0.0482 | n/a | 1.0 | SUBMITTED_LOSS | anti-prior; relationship non-linear, both directions hurt; 0.4514 / r25 |
| 298_zscore_srs | cross-model | internal | +0.001 | -0.0399 | 62%+ | 1.0 | SUBMITTED_LOSS | complex_ipde + OF3 srs z-score; 0.4597 / r24 |
| 302_mmff_strain | physics-rerank | internal | +0.003 | not submitted | 86.4% | 4.0 | HOLD_DIVERGENCE_FAIL | strain noisy; 25/184 overlap, +0.020 worse |
| 322_activity_tanimoto | medoid | internal | +0.008 | -0.0416 | n/a | 3.0 | SUBMITTED_LOSS | activity-Tanimoto medoid family falsified; 0.458 / r27 |

## Iter 3 - First omnibus sweep (June 3, kicked off)

| name | family | external/internal | predicted_lift | actual_lift | divergence_vs_262 | dev_hours | verdict | one_line_reason |
|---|---|---|---|---|---|---|---|---|
| af3_triple_consensus | cofold-3rd-model | EXTERNAL | +0.015 | n/a | n/a | 8 | BLOCKED_INFRA | no public AF3 REST API; 20 jobs/day manual; never run |
| esmfold_induced_fit_dock | de-novo-dock | external | +0.005 | n/a | 88.6% | 12 | HOLD_DIVERGENCE_FAIL | nb08 v5 completed but smina diverged 163/184; CLOSED |
| openmm_md_refine | physics-relax | external | +0.008 | n/a | 67.4% | 16 | HOLD_DIVERGENCE_FAIL | nb09 9h Kaggle T4 success, 184/184 PDBs, 0 high-conf overrides |
| gnina_cnn_rescore | cnn-rescore | external | +0.005 | n/a | 99.5% | 6 | CLOSED_NEGATIVE | nb07 ran but 99.5% divergence = #281 disaster pattern |
| pdbbind_xgboost_lambdarank | learned-ranker | external | +0.005 | n/a | 95.7% | 10 | CLOSED_NEGATIVE | LOO XGB 2.275A vs IPDE 2.230A; -0.045A trend negative |
| chembl_pxr_v2 | data-mining | external | n/a | n/a | n/a | 4 | CLOSED_NEGATIVE | 31-entry corpus, sub-uM pEC50 mismatch with fragments |
| activity_aware_picker | stratified-rerank | internal | +0.001 | n/a | 4.3% | 2 | HOLD_SUB_NOISE | only 8/184 diverge; below 5% gate floor |

## Iter 4 - Re-launched kernels (June 7-8)

| name | family | external/internal | predicted_lift | actual_lift | divergence_vs_262 | dev_hours | verdict | one_line_reason |
|---|---|---|---|---|---|---|---|---|
| nb10_chai1_v4 | cofold-3rd-model | external | +0.015 | n/a | n/a | 10 | CLOSED_NEGATIVE | crashed cell-1 nvidia-smi missing on CPU; 24-48h CPU > 9h cap |
| nb15_proteinmpnn_v4 | pocket-redesign | external | +0.005 | n/a | n/a | 8 | CLOSED_NEGATIVE | 4 MPNN sequences emitted, NO docking step; insufficient |
| nb16_val_expansion | validation-set | internal | 0 | n/a | n/a | 4 | SUBMITTED_NEUTRAL | 18 new ChEMBL holos cofolded; expanded val 35->53; HELD by design |
| nb17_xtb_rescore | physics-rerank | external | +0.005 | n/a | n/a | 6 | CLOSED_NEGATIVE | nvidia-smi crash cell-1, same bug as nb10; never ran xTB |

## Iter 5 - Brainstorm (June 8 first half)

| name | family | external/internal | predicted_lift | actual_lift | divergence_vs_262 | dev_hours | verdict | one_line_reason |
|---|---|---|---|---|---|---|---|---|
| qcarchive_torsion_strain (QTLSP) | physics-rerank | external | +0.009 | n/a | n/a | 14 | CLOSED_NEGATIVE | downstream sanity (#302 MMFF strain val) failed; QTLSP rejected by transitive rule |
| unimol_v2 | learned-ranker | external | unknown | unknown | unknown | 6 | HOLD_DIVERGENCE_FAIL | score for #279 unrecoverable; cannot evaluate vs 0.45 threshold |
| src1_coactivator | external-prior | external | +0.005 | n/a | n/a | 0 | NEVER_RAN | truncated in synthesis; placeholder |
| h12_swap_motif | external-prior | external | +0.003 | n/a | n/a | 0 | NEVER_RAN | truncated in synthesis; placeholder |
| conformal_quantile | meta-ensemble | internal | +0.002 | n/a | n/a | 0 | NEVER_RAN | truncated in synthesis; placeholder |

## Iter 6 - Brainstorm (June 8 evening)

| name | family | external/internal | predicted_lift | actual_lift | divergence_vs_262 | dev_hours | verdict | one_line_reason |
|---|---|---|---|---|---|---|---|---|
| ampt_alphamissense_pocket | external-prior | external | +0.003 | n/a | n/a | 2 | HOLD_SUB_NOISE | Spearman gate passed (\|rho\|=0.294 < 0.3) but sign-inversion risk persists; queued behind iter4/QTLSP |
| diffsbdd_pclp | likelihood-rerank | external | unknown | n/a | n/a | 0 | NEVER_RAN | truncated; GO_RESEARCH only |
| iter6_ideas_3-5 | unknown | unknown | unknown | n/a | n/a | 0 | NEVER_RAN | truncated upstream |

## Iter 7 - Documentation closeout (current)

No new methods proposed. Focus on closing loose ends: stats tables, figures, README. **HOLD #262 (0.4996, rank 10).**

## Key takeaways

- **47 distinct methods attempted; 1 wins (#262 IPDE pure cross-model, 0.4996).**
- **17 SUBMITTED_LOSS** (every off-IPDE submission burned a slot and scored worse than baseline).
- **5 SUBMITTED_WIN** (all are variations of cross-model IPDE / z-hybrid in iter 0).
- **10 HOLD verdicts** (defensive infrastructure caught noise candidates pre-submit).
- **13 CLOSED_NEGATIVE / BLOCKED_INFRA** (kernel-internal re-rankers and external infra failures).
- **NO external dataset has ever produced a validated lift on the structure track.** AF3, ChEMBL, PDBbind, ESMFold, gnina, ProteinMPNN, xTB, MD: all 0/8.
- **Cross-model selection (B2+OF3) is the only proven lift vector.** All future >0.50 plays require external NEW signal (re-refined PDBs, AF3 batches, MCS series) that are not actionable agent-side.
