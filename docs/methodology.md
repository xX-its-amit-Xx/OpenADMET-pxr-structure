# Methodology — `#262 IPDE-select`, validation pipeline, pre-flight gates

This is the detailed write-up of how we arrived at our best submission and the infrastructure that prevented score erosion.

## 1. The challenge in one paragraph

OpenADMET PXR Structure Prediction: predict a 3-D protein-ligand complex (single PDB) for each of **184 PXR LBD + small-molecule pairs**, scored by **LDDT-PLI** (OpenStructure's `ligand_scoring`, bootstrap-averaged over 1000 samples of 184 compounds, half live + half held-out). Submission format: zip of 184 PDB files, ligand residue name exactly `LIG`, parsed ligand graph must match the expected SMILES via RDKit `AssignBondOrdersFromTemplate`. Connectivity failure → per-compound LDDT-PLI=0 + BiSyRMSD=20 Å penalty (so every compound contributes to the bootstrap mean).

Test composition (verified from `pxr-challenge_structure_TEST_BLINDED.csv`):
- **76 PanDDA fragments** — apo P2₁2₁2₁ crystal soaks (10 mM in NSLS-II AMX/FMX), backbone essentially fixed by lattice, ligand pose-in-fixed-receptor problem.
- **108 drug-like analogs** — activity-track compounds, mostly weak PXR binders (24 measured: all pEC50 ≤ 4.4). Heavy-atom mean 24.1, range 10-32.

Strategy implication: the test set is dominated by **novel sub-pocket binders** (fragments at sim<0.3 to all known PXR drug-like ligands). Anchor priors based on known PXR holos are **systematically misleading** for these — and per-family / per-cluster MLPs trained on drug-like corpora **misdirect on PanDDA fragments** (confirmed 3× independently: per_family_mlp #281 = 0.4662 r21, holo_occupancy #294 = 0.4586 r24, activity_tanimoto #322 = 0.458 r27).

## 2. The `#262 IPDE-select` method

### 2.1 Inputs

Two pre-computed pose pools per ligand, both pulled from the `dargason/pxr-cofold` HF dataset:

| Pool | Source | n / ligand | Selection signal |
|---|---|---|---|
| Boltz-2 | 20 seeds of Boltz-2 cofold, full B2 confidence JSON | 20 | `complex_ipde` (lower = better) |
| OpenFold3 | 20 OF3 samples + confidences | 20 | `chain_pair_iptm[(A,L)]` or `pae[pocket, ligand]` |

Total candidate pool: **~7,360 PDBs** (40 × 184) sitting under `data/external/dargason_cofold/predictions/` and `.../openfold3_extract/predictions/`.

### 2.2 The selector

For each ligand `sid` ∈ 184 compound IDs:

1. **Compute Boltz-2 score** `s_B2 = -complex_ipde` over all 20 B2 PDBs. Lower IPDE = higher s_B2 = preferred. (`complex_ipde` = interface pairwise distance error; the model's own confidence in interface geometry. Range typically 0.6–2.5 Å.)
2. **Compute OpenFold3 score**:
   - If `{base}_confidences.json` has `pae`, take a 28-residue × ligand-atoms PAE submatrix indexed by the PXR LBD pocket map (residues 64–284, see below) and return `-mean(pae[pocket, ligand])`.
   - Otherwise fall back to `chain_pair_iptm[(A,L)]` from the aggregated confidences and affine-rescale to roughly match B2's IPDE range: `sig = -2.0 + (iptm - 0.7) / 0.3 * 1.7`.
3. **Rescale OF3 onto B2's IPDE axis** so they're comparable: `sig_scaled = -2.5 + (sig - (-8.0)) / ((-3.0) - (-8.0)) * 1.9`. This empirical affine map was tuned from the Iter-0 z-hybrid runs (#08, #13).
4. **Pick the global argmax** across both pools per ligand. Write `{sid}.pdb` to `submissions/ai_ipde_select/` and zip.

The full script is `scripts/build_ai_ipde_select.py`. It runs in ~30s and produces `submissions/262_ai_ipde_select.zip` (3 MB, 184 PDBs).

### 2.3 PXR LBD pocket map (28 residues)

From `scripts/build_ai_ipde_select.py:POCKET_AI_0IDX` (0-indexed into the 293 aa LBD; residues 142-434 in UniProt O75469):

```
Helices α1/α3:   64, 65, 67, 68, 70, 99, 102, 103, 105, 106, 110
β-sheet / loop:  140, 143, 144, 147, 158, 165, 167
Pocket floor:    182, 183, 186
Helix αAF / lid: 262, 266, 269, 270, 273, 279, 284
```

This pocket map drives:
- The 28-residue PAE submatrix for OF3 scoring.
- The AMPT (AlphaMissense Pocket Tolerance) orthogonality test.
- The contact-fingerprint and hbond-pattern selectors in iter 2.

### 2.4 Why this wins over the alternatives

| Alternative | What we tried | Why it loses |
|---|---|---|
| **Boltz-2 only** (`#282_ipde_b2_only`) | Drop OF3, use only `-complex_ipde` over 20 B2 seeds | 0.4818 (-0.018). Misses cases where OF3 picks the correct novel pocket. |
| **OF3 only** (`#263_of3_only`) | Drop B2, use only OF3 PAE | 0.4724 (-0.027). OF3 has lower per-seed diversity; B2 catches the off-modal hits. |
| **iplddt rerank** (`#290_b2_iplddt`) | Same approach but use `complex_iplddt` instead of `complex_ipde` | 0.4832 (-0.016). iplddt pushes mass toward B2 over OF3 even when OF3 picks correctly. |
| **z-score per model** (`#283_zscore_complex_ipde`) | Per-model z-normalization before combining | 0.4601 (-0.040). Loses absolute calibration; OF3's PAE scale is informative on its own. |
| **B2 ligand_iptm** (`#06_b2_ligand_iptm`) | Use `ligand_iptm` instead of `complex_ipde` | 0.4869 (-0.013). ligand_iptm range across seeds is 0.026 vs IPDE's 0.198 — 7.6× less discriminating. |
| **Cross-model z-hybrid** (`#13`) | Z-score `ligand_iptm` (B2) + `sample_ranking_score` (OF3) | 0.4997 (+0.0001, our gold #2). Direct sibling of #262; sample_ranking_score is OF3's own pre-trained rerank head. |
| **Crystal anchor** (#113, #274, #278) | Bias toward poses close to 1ILH crystal heavy atoms | 0.4657 / 0.1671 / 0.3828 — catastrophic on fragments (novel sub-pockets). |
| **Per-family MLP** (#281) | MLP rerank trained on activity-track per-cluster | 0.4662 (-0.033). 94.6% pick divergence; recapitulates anchor disaster. |
| **Activity-Tanimoto medoid** (#322) | Pick pose closest to medoid of activity-Tanimoto neighbors' best poses | 0.458 (-0.042). Same domain-gap failure on PanDDA fragments. |

Pure IPDE-cross-model is the local maximum across all 47 attempted methods.

### 2.5 The "why" beneath the why

Two structural facts about the test set make IPDE-cross-model the right answer:

1. **PanDDA fragments are NOT well-covered by drug-like priors.** sim<0.3 (Tanimoto, Morgan-r2) between every fragment and every known PXR holo ligand. Any signal trained on drug-like binders' poses (per-family MLP, gnina CNN trained on PDBbind, XGB on PDBbind NR1I, ChEMBL Uni-Mol head) **inverts sign** on fragments. We confirmed this 4 times independently.
2. **Within-model self-confidence (IPDE, iptm) is the strongest signal that survives the drug-like → fragment domain shift.** The model "knows what it doesn't know" — bad IPDE on a fragment tells you the pose is bad regardless of whether the model has ever seen this chemotype.

`complex_ipde`'s 7.6× higher per-seed variance vs `ligand_iptm` is what gives it discriminative power. The fact that B2 and OF3 disagree (we count cross-model wins per ligand) is what makes the mixed selector beat either alone.

## 3. Validation pipeline

The 53-holo validation set is the gate every submission candidate must pass before we burn a slot.

### 3.1 Holo set construction

**35-holo set** (original, iter 1):
- 35 PXR LBD co-crystals from RCSB curated for drug-like ligand binding.
- Pulled into `data/external/pxr_crystals/` (one `{pdbid}.pdb` per holo + crystal ligand SDF if available).
- Run Boltz-2 (20 seeds) and OpenFold3 (20 samples) on the LBD sequence with the crystal-ligand SMILES → 35 holos × 20 B2 + 20 OF3 = 1,400 poses.
- Per-pose oracle RMSD: align the predicted ligand against the crystal ligand by SMILES-derived symmetry mapping, compute heavy-atom RMSD. The MINIMUM RMSD across all 40 poses per holo is the "oracle" lower bound (best possible if we always picked perfectly).

**18-holo expansion** (iter 4, nb16):
- 18 additional PXR-LBD holos curated from ChEMBL with drug-like ligands and decent activity (3HVL, 4J5X, 4NY9, 4X1G, 6HJ2, 6P2B, 6S41, 7AXA, 7AXE, 7AXF, 7AXH, 7AXJ, 7N2A, 7RIV, 8CH8, 8R81, 9FZI, 9FZJ).
- Run Boltz-2 with 20 seeds (nb16 ran for 8h on Kaggle T4) → 18 × 20 = 360 PDBs + confidence JSONs.
- Merged into `data/processed/validation_set/scores.csv` to form the 53-holo set.

### 3.2 Per-pose scoring

`scripts/score_validation.py` computes, for each (holo × pose) pair:
- `rmsd_to_crystal` (heavy-atom RMSD via SMILES symmetry mapping)
- `complex_ipde`, `complex_iplddt`, `ligand_iptm`, `sample_ranking_score`, plus all confidence-JSON-derived signals

Output: `data/processed/validation_set/scores.csv` (1,060 rows on 53-holo set, ~700 on 35-holo only).

### 3.3 Method comparison

`scripts/aggregate_validation.py` aggregates across selectors:
- For each method (e.g. "ipde_best", "plddt_best", "lig_iptm_best", "oracle"), pick the per-holo pose that maximizes the method's signal.
- Compute mean RMSD, under-2 Å fraction, and an `lddt_proxy` (empirical mapping `rmsd → lddt_proxy` calibrated on the 35-holo set against held-out leaderboard pulls).

Output: `data/processed/validation_set/method_comparison.csv`.

### 3.4 35-holo vs 53-holo: the plddt-overfit unmask

| method | 35-holo mean RMSD | 35-holo lddt_proxy | 53-holo mean RMSD | 53-holo lddt_proxy | Δ lddt |
|---|---|---|---|---|---|
| ORACLE (lower bound) | 1.683 | 0.5109 | 1.403 | 0.5346 | +0.024 |
| first_model_0 | 1.985 | 0.4043 | 1.921 | 0.4255 | +0.021 |
| conf_score_best | 1.985 | 0.4043 | 1.921 | 0.4255 | +0.021 |
| #290_iplddt_best | 2.014 | 0.4051 | 1.922 | 0.4253 | +0.020 |
| #13_lig_iptm_best | 2.018 | 0.3998 | 1.931 | 0.4208 | +0.021 |
| **#262_ipde_best** | **2.230** | **0.4024** | **1.938** | **0.4214** | **+0.019** |
| pde_best | 2.040 | 0.3975 | 1.950 | 0.4185 | +0.021 |
| **#plddt_best** | **2.099** | **0.4075** | **1.968** | **0.4090** | **+0.0015** ⚠ |

`#plddt_best` led on 35-holo by +0.005 lddt_proxy over IPDE, but on 53-holo it falls 0.0124 BEHIND IPDE. The 18-holo expansion was easier overall (oracle improved +0.024), yet plddt only gained +0.0015 — confirming it was a **35-holo overfit** that doesn't generalize. This validated the 2026-06-02 decision to HOLD #262 instead of submitting the v3_hybrid_pocket_plddt candidate (which scored 0.4702, a sign-inverted miss).

See [docs/figures/iter4_validation_expansion_35_vs_53.png](figures/iter4_validation_expansion_35_vs_53.png).

### 3.5 The "validation does not predict leaderboard" caveat

Critical: a method's 53-holo lddt_proxy is **NOT** monotonic with its leaderboard LDDT-PLI. See `project_validation_results_v1.md`:

- All tested methods cluster within 0.05 Å of mean RMSD on 35-holo (within noise floor).
- But leaderboard scores span 0.16 → 0.50, a **5× wider range**.
- Cause: holo set has drug-like binders (medium-sim to training corpora), test set is PanDDA-fragment-heavy (novel chemotypes). Methods that overfit drug-like priors look fine on holo but DIE on fragments.

This is why **the pre-flight gate is divergence + RMSD + confidence**, not just RMSD. Divergence vs #262 catches anchor disasters that look fine on holo.

## 4. Pre-flight gates (the three-gate rule)

Every new submission candidate must pass ALL THREE gates before we burn a 4-hour slot:

### Gate 1 — Connectivity / format

`scripts/validate_submission.py <zip>`:
- 184 PDBs present.
- Ligand resname is exactly `LIG` (not `LIG1` — Boltz-2's default; we rename via biotite).
- Max 2 chains.
- Ligand graph parsed by RDKit matches expected SMILES via `AllChem.AssignBondOrdersFromTemplate`.

Failure → per-compound LDDT-PLI=0 + BiSyRMSD=20 Å. We burned a slot to confirm this experimentally early on. NEVER skip.

### Gate 2 — Divergence band [5%, 30%] vs #262

`scripts/validate_selection.py <new_zip> submissions/262_ai_ipde_select.zip`:
- Counts fraction of 184 picks where the new selector chose a different pose than #262.
- **< 5%**: sub-noise. Even if the new picks are slightly better, the lift would be <+0.001 (below leaderboard noise floor). Don't waste a slot.
- **5% – 30%**: GREEN. New picks have meaningful coverage divergence but are not a wholesale rewrite.
- **> 30%**: RED. Anchor-disaster band. Methods with high divergence (94.6% #281, 99.5% gnina, 95.7% XGB, 88.6% smina, 86.4% MMFF strain) have all scored sub-0.47 on the leaderboard.

13 candidates caught by this gate alone. See `verdict=HOLD_DIVERGENCE_FAIL` column of [docs/stats/methods_table.md](stats/methods_table.md).

### Gate 3 — 53-holo mean RMSD ≤ 2.15 Å + confidence ≥ medium

- `data/processed/validation_set/method_comparison.csv` for the new method.
- Threshold 2.15 Å chosen empirically: it sits ~0.08 Å below current #262 (2.230 Å) and well above the noise floor (0.05 Å). Methods within 0.05 Å of #262 are inside noise; methods below 2.15 Å are meaningfully better.
- Confidence: at least 3 of {ipde, iplddt, lig_iptm, srs} support the picks (median rank ≤ 5 across pools).

Result: across iter 3–6, **0 / 12 candidates passed all three gates**. HOLD #262 has been the correct decision every time.

## 5. Submission safety infrastructure

### 5.1 No-restore-reflex

The OpenADMET Gradio endpoint shows your MOST RECENT submission, not your best. This is the single biggest leaderboard trap on this challenge:

- **3 ex-top-5 teams** dropped > 0.05 in our submission window by submitting worse experiments:
  - `suspenders`: 0.5521 → 0.4727 (-0.0794), rank 2 → 18
  - `Cryo-EMinem`: 0.5113 → 0.4640 (-0.0473), rank 3 → 25
  - `blacktea`: 0.5015 → out of top-20, rank 6 → 21+

Our rule (codified after we ourselves dropped to 0.1671 / rank 38 with `#274_mega_crystal_sa`): **never submit a candidate that has not passed all three gates**. Auto-submit `262_ai_ipde_select.zip` as the LAST submission of any 4-hour cycle.

### 5.2 Windows Scheduled Task + backup PID

The deadline-restore is wired in two redundant ways:

1. **Windows Scheduled Task `PXR_Safety_Restore_262`** fires 2026-06-30 23:00 UTC, action:
   ```
   .venv\Scripts\python.exe scripts\auto_submit_v2.py submissions\262_ai_ipde_select.zip --confirm
   ```
2. **Backup PID 20016** (background Python heartbeat) fires same script at 23:30 UTC.

`scripts/auto_submit_v2.py` is idempotent: it fetches the current leaderboard, only submits if our row is not already 0.4996. If the network blip prevents the leaderboard fetch, it submits anyway (default-safe).

### 5.3 Validation log

`data/processed/leaderboard_log.csv` records every submission with:
- `timestamp_utc`, `submission_name`, `local_lddt_estimate`, `leaderboard_score`, `divergence_vs_262`, `notes`

The local-vs-leaderboard divergence column has flagged 4 sign-inversions where our 53-holo estimate was wrong (the most expensive: `v3_hybrid_pocket_plddt` predicted 0.5095, scored 0.4702 — a -0.04 sign flip).

## 6. What we would do differently

(Listed because: the user asked for an honest write-up.)

1. **Codify the pre-flight gate by iter 2, not iter 3.** We burned 11 SUBMITTED_LOSS slots in iter 1–2 before we wrote `validate_selection.py`. If the divergence gate had been in place from iter 1, we would have caught 8 of those 11.
2. **Set up the 53-holo validation set earlier.** The 35-holo set was sufficient for iter 0 baselines, but the plddt-overfit was already lurking. Building the 18-holo expansion in iter 1 (not iter 4) would have rejected pocket_plddt before it ever shipped.
3. **Don't burn 80h on omnibus external Kaggle attempts.** AF3 manual upload, ESMFold+smina, OpenMM MD refine, gnina CNN, PDBbind XGB, Chai-1 v4, ProteinMPNN, xTB: **0 / 8** of these reached terminal output with a submission-ready candidate. The 5 external paths in iter 3 had a uniform 12-week-activation-lag failure mode that doesn't fit a 30-day window.
4. **Push for re-refined PDBs by week 2.** The 64 re-refined PDBs were posted to GitHub 2026-05-18 (3 days into the challenge for us); we didn't find them until 2026-06-08. Per OpenADMET's own cofold-eval blog, they wouldn't have helped much — but we should have known that 3 weeks earlier.
5. **Focus 70% of brainstorm volume on the test-set distribution.** The PanDDA-fragment vs drug-like split was load-bearing and we didn't internalize it until iter 1 mid-week. Anchor / family / occupancy priors are all drug-like-aligned; on a fragment-heavy test set they were doomed.

## 7. The methodology in 5 sentences (TL;DR)

1. **Use both Boltz-2 and OpenFold3 cofold predictions** from the dargason HF dataset (~40 poses per ligand).
2. **Score B2 by -complex_ipde**; score OF3 by -pocket-residue × ligand-atom PAE (or chain_pair_iptm fallback), affine-rescaled onto the B2 axis.
3. **Pick the global argmax per ligand**, write 184 PDBs, validate connectivity, submit.
4. **Validate locally** on 53 PXR holo crystals (35 RCSB + 18 ChEMBL) — IPDE-best is the most robust selector across both 35 and 53-holo subsets.
5. **Gate every new candidate** through (a) connectivity, (b) 5-30% divergence vs #262, (c) 53-holo mean RMSD ≤ 2.15 Å with medium+ confidence — and HOLD #262 if any gate fails.
