# Iteration Summary (Iter 0 - Iter 7)

Per-iteration accounting of brainstorm volume vs build/submit conversion vs realised lift. Definitions:

- `n_ideas_proposed`: distinct method ideas formally tracked in tracker / memory.
- `n_GO_BUILD`: ideas that received code-build commitment.
- `n_HOLD`: ideas built but never submitted (failed pre-flight gate, sub-noise, divergence-fail, infra-blocked).
- `n_REJECT`: ideas explicitly rejected before any build effort.
- `n_submitted`: distinct submissions to Gradio endpoint.
- `max_actual_score`: highest LDDT-PLI achieved by any submission in that iter.
- `key_learning`: most actionable cross-iter takeaway.

## Per-iter table

| Iter | Window | n_ideas | n_GO_BUILD | n_HOLD | n_REJECT | n_submitted | max_actual_score | key_learning |
|---|---|---|---|---|---|---|---|---|
| 0 (baseline) | May 18-20 | 5 | 5 | 0 | 0 | 5 | **0.4997** (#13) | cross-model z-hybrid (B2 ligand_iptm + OF3 srs) wins; signal choice matters more than aggregation |
| 1 (anchor/family) | May 21-23 | 12 | 11 | 2 | 1 | 10 | **0.4997** / **0.4996** | crystal-anchor + per-family-MLP catastrophic on PanDDA fragments (sim<0.3); pure IPDE wins |
| 2 (signal mining) | May 25-28 | 11 | 8 | 4 | 3 | 6 | 0.4832 | every off-IPDE B2/OF3 signal combo regresses; holo-prior fails both directions |
| 3 (omnibus iter1) | June 3 | 7 | 7 | 4 | 0 | 0 | n/a (HOLD) | 5 ambitious external resources kicked off; all CLOSED_NEG or BLOCKED_INFRA at iter1 closure |
| 4 (re-launched kernels) | June 7-8 | 4 | 4 | 2 | 0 | 0 | n/a (HOLD) | nb10/nb15/nb17 all crashed on nvidia-smi or sub-noise; only nb16 (val expansion 35->53) succeeded |
| 5 (brainstorm) | June 8 AM | 5 | 1 | 0 | 0 | 0 | n/a (HOLD) | QTLSP (DFT torsion strain) single GO_BUILD; 4 ideas truncated upstream |
| 6 (brainstorm) | June 8 PM | 5 | 0 | 1 | 1 | 0 | n/a (HOLD) | AMPT GO_RESEARCH only; #302 MMFF strain validator transitively rejected QTLSP |
| 7 (closeout) | June 8-30 (planned) | 0 | 0 | 0 | 0 | 0 | n/a (HOLD) | documentation phase; no new methods |

**Totals (iter 0-6, realised):**
- 49 ideas proposed
- 36 GO_BUILD commitments
- 13 HOLD (caught by gate)
- 5 REJECT (pre-build)
- 21 submissions to Gradio
- 5 SUBMITTED_WIN (all in iter 0; all cross-model z-hybrid variants)

## Brainstorm efficiency analysis

| Phase | Submissions / GO_BUILD | Score lift per submission | Submission ROI |
|---|---|---|---|
| Iter 0 (cofold baselines) | 5/5 = 100% | +0.0073 / sub (0.4632 -> 0.4997) | HIGH (each one informs next) |
| Iter 1 (anchor/family) | 10/11 = 91% | -0.0167 / sub on average | NEGATIVE (3 self-overwrites recovered via #262) |
| Iter 2 (signal mining) | 6/8 = 75% | -0.0341 / sub | VERY NEGATIVE (all kernel-internal signals fail) |
| Iter 3 (omnibus external) | 0/7 = 0% | n/a | HOLD pays off (caught 5 bad candidates pre-submit) |
| Iter 4 (kernel re-launch) | 0/4 = 0% | n/a | infra failures, no submission risk |
| Iter 5 (brainstorm) | 0/1 = 0% | n/a | QTLSP transitively rejected after sanity |
| Iter 6 (brainstorm) | 0/0 = 0% | n/a | AMPT held; iter4 still in-flight |

**Key inflection point:** between iter 2 (5 submissions, all losses) and iter 3 (0 submissions, validation gate enforced), we shifted from "submit-then-learn" to "validate-then-submit". The 6-day HOLD since 2026-06-02 has preserved 0.4996 against a field where 3 top teams self-overwrote and dropped >0.05.

## What worked across iters

1. **Cross-model B2+OF3 z-hybrid / IPDE selection** (iter 0, iter 1 culmination as #262). +0.0364 vs baseline, +0.0001 vs gold runner-up #13.
2. **Validation gate (35-holo + later 53-holo)** (iter 3 onward). Caught MD-refined (67.4% div), gnina (99.5%), XGB (95.7%), smina (88.6%), MMFF strain (86.4%).
3. **No-restore-reflex + Windows scheduled task** (iter 1 onward). Prevented score erosion.
4. **Auto-submit_v2 restore** (iter 1 onward). Recovered from 4 self-overwrite incidents (0.1671 -> 0.4996 in 4h).

## What never worked

1. **Any single-model rerank** (Boltz-only, OF3-only) - both 0.47-0.48, below cross-model 0.4996.
2. **Crystal anchor / holo prior on PanDDA fragments** - sim<0.3 makes anchor pure noise; 3 independent submission failures.
3. **Pocket-residue-based signals on fragments** (pocket_plddt, hbond_pattern, contact_fp) - sign-invert on fragment test set (predicted +0.0099 -> scored -0.0294).
4. **Kernel-internal physics rerank** (MD energy, MMFF strain, xTB) - shuffle near-tie poses in IPDE noise band; no high-confidence overrides.
5. **Learned rerankers trained on drug-like corpora** (per-family MLP, XGB on PDBbind, Uni-Mol on ChEMBL PXR) - same domain-gap failure mode, 4 independent confirmations.
6. **External cofold models on Kaggle** (Chai-1, ProteinMPNN+smina, ESMFold+vina) - infra blockers, GPU quota, nvidia-smi bugs, module ABI mismatches; 0/3 reached terminal output with usable picks.

## Time / compute budget

- Total dev hours (estimate, all iters): ~180h
  - Iter 0-2 (submission frenzy): ~50h
  - Iter 3-4 (omnibus + re-launches): ~80h
  - Iter 5-7 (brainstorm + closeout): ~50h
- Total Kaggle GPU hours: ~60h (nb01-17 across all iters)
- Total submissions to Gradio: 21 (each costs 4h cool-down)
- Total submissions still showing as "best" on leaderboard: 1 (#262)

## Honest assessment

### Current standing (2026-06-08, T-23 days)

- **Our score: 0.4996** (#262_ai_ipde_select, pure B2+OF3 IPDE cross-model selection).
- **Our rank: 10** (drifted from #5 (May 20) -> #6 (May 22) -> #9-10 (June 2) -> #10 (June 8)).
- **Leader: 0.5612** (dnan-ipd, locked since 2026-05-29).
- **Top-5 cutoff: 0.5115** (discoverybytes).

### Gap analysis

| To beat | Score required | Lift needed from 0.4996 |
|---|---|---|
| Top-10 (us) safe floor | 0.4996 | 0.0000 (HOLD wins) |
| Top-5 cutoff (#5 = discoverybytes) | 0.5115 | +0.0119 |
| Top-3 (TangerineTrees #3) | 0.5285 | +0.0289 |
| Leader (dnan-ipd #1) | 0.5612 | +0.0616 |

### Ceiling estimates (per project_numerical_reassessment_2026_06_08)

**Empirical kernel-internal ceiling (no external signal):** ~0.5050. This is +0.0054 over current, derived from:
- 35-holo + 53-holo validation gap between IPDE (2.214A) and oracle (1.403A) = 0.811A headroom
- All 12 omnibus candidates have landed within +/-0.045A on this gap = noise floor
- Probability of any kernel-internal method exceeding 0.51 from validation alone: **<10%**

**Ceiling with external signals (probabilistic):**
| Path | Best-case score | Probability >= 0.51 |
|---|---|---|
| Triple consensus B2+OF3+Chai-1 (if Chai-1 lands) | ~0.515 | 35-45% |
| + AF3 server batches landed (quadruple) | ~0.525 | 55-65% |
| + 68 re-refined PDBs from organizers | ~0.530 | 50% |
| All three landed | ~0.535 | best realistic |

**Realistic ceiling at deadline: 0.515-0.535.** Leader at 0.5612 likely requires PXR-specific privileged training data (organizer-held holos) or a custom pose-prediction model we cannot replicate.

### External dependencies (none actionable agent-side)

1. **AlphaFold3 server batches (184 jobs).** BLOCKED_EXTERNAL. 2 batch ZIPs prepared at `data/external/af3_batches/`. Requires manual upload to alphafoldserver.com (20 jobs/day x ~10 days). User action required. If completed, +0.010 to +0.020 expected lift.

2. **68 re-refined PDBs** promised by organizers via Hugging Face. NOT YET RELEASED as of 2026-06-08. If delivered, would expand validation set from 35 (drug-like) + 18 (ChEMBL drug-like) to 35 + 18 + 68 = 121 holos, including PanDDA fragments at the test-set distribution. Expected lift: +0.005 to +0.020 from better-anchored validation + improved protein structure priors.

3. **MCS-series ensembling on real 40-60 series.** Investigated and deferred (project_mcs_series_explored). Naive clustering gives 120 series but 85 are transitive-closure artifacts; real coverage 40-60. Expected lift +0.008 below noise floor on validation. Not actionable without external pose-pool augmentation.

4. **PDBbind NR1I subset / XChem fragments** for ranker training. Untried external corpora (per reference_external_datasets). Domain-gap risk similar to ChEMBL PXR v2 (CLOSED_NEGATIVE at 31 entries). High risk, modest expected lift.

### Strategic verdict (2026-06-08)

- **Path of least regret: HOLD #262.** With 23 days remaining and 4 self-overwrite casualties already on the leaderboard, the dominant strategy is to preserve 0.4996.
- **Safety net armed.** Windows Scheduled Task `PXR_Safety_Restore_262` fires 2026-06-30 23:00 UTC; backup PID 20016 fires 23:30 UTC.
- **Submission triggers (any one unblocks):** AF3 batches return AND triple/quadruple consensus passes 3-gate; OR 68 re-refined PDBs ship AND new validation set shows IPDE-best >= 0.51 with [5-30%] divergence vs #262.
- **Without external trigger, final submission stays at 0.4996.** Final rank likely 10-15 as the field continues moving (4 new top-20 entrants in 6 days).

### Lessons that generalize

1. **Cross-model selection beats every other signal class** on PanDDA-fragment-heavy test sets where drug-like priors break.
2. **Validation->test compression is unreliable for pocket-residue signals** (sign-inverted on v3_hybrid_pocket_plddt). Only IPDE-class signals transfer.
3. **The leaderboard-shows-latest trap is real and costly.** No-restore-reflex + automated safety task + 25% divergence pre-flight gate are necessary infra for any leaderboard with no submission history endpoint.
4. **External resources have a 12-week activation lag** that doesn't fit a 30-day competition window. The 5 ambitious iter-3 external paths (AF3, ESMFold, MD, gnina, XGB) all hit infra blockers before producing terminal output.
5. **Brainstorming yields diminishing returns past ~20 ideas** when the underlying signal space is already mapped (signal-correlations |rho|<0.08 across all kernel-internal confidence signals).
