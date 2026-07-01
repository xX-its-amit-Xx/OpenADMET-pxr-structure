# OpenADMET PXR — Structure Prediction Track

Final documentation for our entry in the **OpenADMET PXR Blind Structure-Prediction Challenge** (2026-05-17 → 2026-07-01).

---

## Final standing

| Metric | Value |
|---|---|
| **Best score** | **0.5640 LDDT-PLI** |
| **Best rank** | **2 / ~50 teams** |
| **Best submission** | `prot_rescue8` |
| Leader | 0.5725 (`dnan-ipd` — Apheris AI, federated PL fine-tune) |
| Gap to leader | +0.0085 |
| Challenge close | 2026-07-01 ~23:59 UTC |

The leaderboard shows the **most recent** submission, not the best. Our OS Scheduled Task (`pxr_auto_submit_structure`) fires at 2026-07-01 22:30 UTC to guarantee `prot_rescue8.zip` is the final displayed submission.

---

## Winning method: `prot_rescue8`

### One-line description

4-model z-hybrid confidence selector (AF3 + Boltz-2 + OpenFold3 + Chai) grafted with Protenix-v2 poses on the 8 lowest-confidence ligands.

### How it works

**Step 1 — Multi-model pose pool**

For each of the 184 test compounds we have a pool of candidate poses from 6 co-folders:

| Model | Source | # seeds/samples |
|---|---|---|
| AlphaFold3 (AF3) | Modal A100 + alphafoldserver.com | 4–10 per ligand |
| Boltz-2.1 (B2) | Boltz hosted API ($100 credits) | 5–10 per ligand |
| OpenFold3 (OF3) | dargason/cofold HF dataset | 20 per ligand |
| Chai-1 (Chai) | Explorer cluster / Colab Pro | 5 per ligand |
| Protenix-v2 | Explorer cluster (SLURM 7762644) | 25 per ligand |
| ESMFold2 | Explorer cluster (4 MSA modes) | 25 per mode |

Pool oracle (best achievable RMSD): ~1.08 Å median across all 184 ligands.

**Step 2 — Within-model z-score selection**

For each model independently we pick the best seed/sample using per-model confidence signals z-scored within that model's pool:

- **AF3:** `iptm` (interface pTM)
- **Boltz-2:** `-complex_ipde` (negative complex interface PDE)
- **OpenFold3:** `-pae_pocket` (negative mean pocket-ligand PAE)
- **Chai:** `iptm`

Each model contributes one "best" pose per ligand after within-model selection.

**Step 3 — Cross-model z-hybrid selection**

Across the 4 base models, we pick the model with the highest z-scored confidence signal for each ligand independently. This is the core breakthrough: cross-model diversity (each model fails on different ligands) combined with z-normalization (prevents one model's inflated raw scores from dominating) produces the best per-ligand pick.

Exact scoring recipe:

1. For every model `m`, ligand `i`, and sample `k`, compute a model-native raw confidence `r[m,i,k]`.
2. Pick the best sample inside that model: `k* = argmax_k r[m,i,k]`.
3. Z-score the best-sample scores across all 184 ligands for that model:
   `z[m,i] = (r[m,i,k*] - mean_i r[m,i,k*]) / std_i r[m,i,k*]`, with `std=1.0` as the zero-variance fallback.
4. Copy the pose from `argmax_m z[m,i]`.

Model signals:

| Model | Signal |
|---|---|
| Boltz-2 | `-complex_ipde` |
| OpenFold3 | `-mean(PAE[protein residues 1:293, ligand tokens])`, fallback `chain_pair_iptm(A,L)` |
| AF3 | exported AF3 ranking/interface score from `af3_best/_ranking.json` |
| Chai-1 | exported Chai ranking/interface score from `chai_best/_ranking.json` |

Script: `scripts/build_ai_4model_zhybrid.py` — base score: **0.5472**

**Step 4 — Protenix-v2 tail rescue**

Protenix-v2's 25-sample pool dramatically rescues the deepest failure-mode ligands (example: 4X1G ligand 0.123 → 0.919, +0.796). We identify the 8 ligands with the lowest cross-model confidence z-score (the "failure tail") and overwrite their z-hybrid picks with the best Protenix-v2 pose for those 8 ligands.

Swap bracket validation (N=number of Protenix swaps):

| N | Score |
|---|---|
| 4 | 0.5578 |
| **8** | **0.5640 ← peak** |
| 12 | 0.5629 |
| 20 | 0.5587 |

8 swapped ligands: `x00035, x00242, x00337, x00558, x00990, x01131, x01334, x01438`

Script: `scripts/build_prot_rescue.py 8 prot_rescue8` — final score: **0.5640**

**Step 5 — Validation + submission hygiene**

Before every submission:
1. `scripts/validate_submission.py` — PDB format, `LIG` resname, ligand-graph SMILES match (zero errors required)
2. `scripts/clean_scoring_fails.py` — swap 3 server-confirmed scoring-fail ligands (`x03063, x03260, x03264`) to clean fallback poses; the full strict check over-flags, so we clean only confirmed failures
3. `scripts/add_conect.py` — inject CONECT records into all 184 PDBs so the server infers ligand bonds from topology, not 3D geometry

---

## Score history

All realized leaderboard scores, ordered by submission date:

| Date (UTC) | Tag | Score | Rank | Note |
|---|---|---|---|---|
| 2026-05-22 | `262_ai_ipde_select` | 0.4996 | 5→10 | Baseline IPDE selector; HOLD for 6 days |
| 2026-06-21 | `prot5` | 0.5241 | 6 | First Protenix test; 5 swaps |
| 2026-06-22 | `sub_af3_chai_b2` | 0.5472 | 3 | 4-model z-hybrid base (AF3+Chai+B2) |
| 2026-06-22 | `sub_af3_of3` | 0.5227 | 7 | 2-model AF3+OF3 |
| 2026-06-22 | `sub_af3_b2` | 0.5405 | 4 | 2-model AF3+B2 |
| 2026-06-22 | `sub_af3_chai` | 0.5212 | 6 | 2-model AF3+Chai |
| 2026-06-22 | `prot_rescue` | 0.5629 | 2 | 12-swap Protenix rescue |
| 2026-06-22 | `prot_rescue20` | 0.5587 | 2 | 20-swap (over-swaps) |
| **2026-06-23** | **`prot_rescue8`** | **0.5640** | **2** | **← OUR BEST** |
| 2026-06-23 | `prot_rescue4` | 0.5578 | 2 | 4-swap |
| 2026-06-23 | `af3deep_upgrade` | 0.5632 | 2 | AF3 deep sampling upgrade |
| 2026-06-23 | `prot_rescue_hs8` | 0.5608 | 2 | Protenix high-sample re-run |
| 2026-06-23 | `af3_iface_upgrade` | 0.5625 | 2 | AF3 interface-signal upgrade |
| 2026-06-23 | `zhyb5_prot` | 0.5241 | 9 | 5-model z-hybrid (regressed) |
| 2026-06-23 | `esmfold2_full` | 0.5264 | 7 | ESMFold2 full-MSA baseline |
| 2026-06-24 | `esmfold2_rescue8` | 0.5610 | 2 | ESMFold2 swap on 8 tail ligands |
| 2026-06-24 | `rmsdpred_rescue8` | 0.5574 | 2 | RMSD-pred learned scorer swap |
| 2026-06-24 | `combo_prot8_esm8` | 0.5563 | 2 | Protenix+ESMFold2 combined |
| 2026-06-26 | `manual_annot` | 0.5558 | 3 | 3-variant medchem annotation |
| 2026-06-26 | `manual_annot_regular` | 0.5613 | 2 | Medchem light/medium swaps only |
| 2026-06-26 | `manual_annot_aggressive` | 0.5578 | 3 | Including drastic swaps (over-corrects) |
| 2026-06-27 | `combo_prot8_boltz1_8` | 0.5362 | 7 | Protenix+Boltz-1 combined |
| 2026-06-27 | `af3_iface_rescue8` | 0.5557 | 3 | AF3 interface rescue |
| 2026-06-28 | `hardtail_rescue` | 0.5456 | 4 | Tail-hardest-10 swap |
| 2026-06-28 | `hardtail20_rescue` | 0.5441 | 5 | Tail-hardest-20 swap |
| 2026-06-28 | `posescorer_pick` | 0.4762 | 32 | Custom XGBoost pose scorer (failed) |
| 2026-06-29 | `threevar_recommended_safe` | 0.5462 | 6 | Medchem annotation safe set |
| 2026-06-29 | `threevar_all_light` | 0.5570 | 4 | All-light medchem edits |
| 2026-06-29 | `flip_moonshot_agentic` | 0.5627 | 2 | Edge-chemistry flip + agentic filter |
| 2026-06-29 | `threevar_recommended_full` | 0.5363 | 10 | Including drastic medchem |
| 2026-06-30 | `agentic_pool_pick` | 0.5140 | 19 | Full agentic PXR-expert selector |
| 2026-06-30 | `prot_rescue8_tail5medchem` | 0.5611 | 3 | Medchem on 5 tail ligands |
| 2026-06-30 | `intfold_base_rescue8` | *(pending)* | — | IntFold base + 8 Protenix swaps |

---

## Repository layout

```
OpenADMET-pxr-structure/
├── README.md                            # this file
├── pyproject.toml / uv.lock             # Python deps (uv-managed)
├── .env                                 # credentials (gitignored)
│
├── scripts/                             # all build/validate/submit/monitor code
│   │
│   │  ── SUBMISSION PIPELINE ──
│   ├── submit.py                        # Gradio /submit_predictions wrapper
│   ├── ladder_submit.py                 # THE single submission authority (priority queue, 4h cooldown)
│   ├── fetch_leaderboard.py             # Live leaderboard fetch + best_submission.json update
│   ├── clean_scoring_fails.py           # Swap 3 confirmed server-scoring-fail ligands before submit
│   ├── add_conect.py                    # Inject CONECT records (bond topology, avoids coord-misparse)
│   ├── apply_overrides.py               # Force-apply hand-built override PDBs before submit
│   ├── validate_submission.py           # PDB format + ligand-graph SMILES match (184-file)
│   ├── validate_selection.py            # Pre-flight divergence gate vs reference submission
│   ├── deadline_restore.py              # OS-task helper: restore prot_rescue8 at deadline
│   │
│   │  ── WINNING METHOD ──
│   ├── build_ai_4model_zhybrid.py       # 4-model z-hybrid base (AF3+B2+OF3+Chai) → 0.5472
│   ├── build_prot_rescue.py             # Protenix tail-swap on N lowest-confidence ligands → 0.5640
│   │
│   │  ── POSE GENERATION ──
│   ├── modal_af3.py                     # AlphaFold3 on Modal A100 (4-seed loop)
│   ├── boltz_api_msa.py                 # Boltz-2.1 via hosted API, full-MSA
│   ├── boltz_api_nomsa.py               # Boltz-2.1 via hosted API, no-MSA
│   ├── boltz_druglike_massive.py        # Best-of-10 Boltz-2.1 on 87 drug-like analogs
│   ├── boltz_gt_validate.py             # Boltz-2.1 best-of-5 on holo crystals (GT gate)
│   ├── build_af3_user_batches.py        # AF3 server batch JSONs for manual upload
│   │
│   │  ── SELECTION & SCORING EXPERIMENTS (all refuted, kept for audit) ──
│   ├── build_ai_3model_zhybrid.py       # 3-model (B2+OF3+AF3) z-hybrid → 0.5414
│   ├── build_ai_5model_zhybrid.py       # 5-model z-hybrid → 0.5241 (regressed)
│   ├── build_adaptive_ensemble.py       # Adaptive signal weighting
│   ├── build_ai_boltzmann_avg.py        # Boltzmann-averaged pose
│   ├── build_ai_consensus_centroid.py   # Multi-model consensus medoid
│   ├── build_adaptive_soft_medoid.py    # Soft-medoid consensus
│   ├── build_ai_contact_oracle.py       # Pocket-contact oracle selector
│   ├── build_ai_conf_pca.py             # PCA of confidence signals
│   ├── build_ai_activity_tanimoto.py    # Activity-Tanimoto medoid selector
│   ├── build_ai_crystal_anchor_sa.py    # Crystal-anchor selection
│   ├── pxr_pose_score.py                # Custom XGBoost LambdaMART pose scorer → 0.4762
│   ├── flip_tester.py                   # Interaction-edge-chemistry flip/rotation search
│   ├── pose_recombine.py                # Genetic anchor+tail pose crossover
│   ├── orient_fix.py                    # Rotational search at fixed centroid
│   ├── local_refine.py                  # OpenMM local MD refinement (refuted)
│   │
│   │  ── VALIDATION & ANALYSIS ──
│   ├── pose_lib.py                      # Reusable GT harness: parse_pose, gt_min_rmsd
│   ├── score_validation.py              # Per-pose RMSD on holo crystal validation set
│   ├── aggregate_validation.py          # Method comparison on 35/53-holo set
│   ├── validate_selection.py            # 5–30% divergence gate
│   ├── ligand_chem_analyze.py           # Per-ligand medchem dossier + RDKit featurizer
│   ├── af3_vs_boltz_gt.py               # AF3 vs B2 head-to-head on holo GT
│   ├── aggregate_validation.py          # method_comparison.csv aggregator
│   │
│   │  ── POST-CHALLENGE (not run during competition) ──
│   ├── finetune/
│   │   ├── prep_sair_nr.py              # Filter SAIR 5.2M cofolds for NR regime
│   │   ├── protenix_lora.sbatch         # Explorer Slurm LoRA fine-tune job (8h chunks)
│   │   ├── finetune_config.yaml         # LoRA config (freeze input_embedder+msa_module)
│   │   ├── eval_ft_on_holos.py          # THE GATE: FT vs public on holo GT
│   │   └── README.md                    # Run order, cost ($0 on Explorer), honest gaps
│   │
│   └── activity/                        # Activity-track scripts (sister repo interface)
│
├── notebooks/                           # Kaggle GPU kernels (historical)
│   ├── 02_boltz2_kaggle_full.ipynb      # Boltz-2 baseline
│   ├── 06_boltz2_validation.ipynb       # 35-holo oracle
│   ├── 16_boltz_chembl_validation.ipynb # 18 ChEMBL holos → 53-holo val set
│   └── … (17 notebooks total)
│
├── submissions/                         # Built submission zips (gitignored)
│   ├── prot_rescue8.zip                 # OUR BEST (0.5640, rank 2)
│   ├── af3_batches/                     # JSONs for manual AF3 server upload
│   └── … (30+ candidate zips)
│
├── data/
│   ├── external/                        # Input datasets (gitignored)
│   │   ├── dargason_cofold/             #   3,680 Boltz-2 + 3,680 OF3 PDBs
│   │   ├── pxr_holo/                    #   9 cached PXR LBD holos
│   │   ├── pxr_crystals/               #   35-holo original validation set
│   │   ├── pxr_crystals_chembl/         #   18 ChEMBL holos (nb16 expansion)
│   │   └── druglike_massive_best/       #   Best Boltz-2 pose per drug-like analog
│   └── processed/
│       ├── ladder_queue.csv             # Priority-ordered submission queue + realized scores
│       ├── ladder_state.json            # Last submit time + in-flight tag
│       ├── best_submission.json         # Canonical best: prot_rescue8 @ 0.5640
│       ├── leaderboard_log.csv          # Every submission attempt + Gradio response
│       ├── validation_set/              # 53-holo RMSD results + method_comparison.csv
│       ├── medchem_kb/                  # PXR SAR/interactions/gotchas knowledge base
│       └── manager_log.md               # Running DID log of all actions this session
│
└── docs/                                # Closeout documentation
    ├── INDEX.md
    ├── methodology.md
    ├── iterations.md
    └── figures/ + stats/
```

---

## What we tried — comprehensive method table

### Successful methods (scored ≥ 0.56)

| Method | Score | Key idea |
|---|---|---|
| `prot_rescue8` | **0.5640** | 4-model z-hybrid + Protenix-v2 on 8 failure-tail ligands |
| `prot_rescue` | 0.5629 | Same but 12 Protenix swaps (slightly over-swaps) |
| `af3deep_upgrade` | 0.5632 | AF3 with deeper sampling per ligand |
| `af3_iface_upgrade` | 0.5625 | AF3 re-scored by interface signal |
| `flip_moonshot_agentic` | 0.5627 | Edge-chemistry flip search, agentic-filtered 6/90 flips |
| `manual_annot_regular` | 0.5613 | 3-variant medchem annotation, light+medium swaps only |

### Methods that regressed (score < 0.564)

| Method | Score | Why it failed |
|---|---|---|
| `prot_rescue_hs8` | 0.5608 | Protenix high-sample rerun; same 8 ligands, marginal |
| `esmfold2_rescue8` | 0.5610 | ESMFold2 on 8 tail; better than 4th-model swap |
| `prot_rescue8_tail5medchem` | 0.5611 | Medchem on 5 tail; medchem edits break clean poses |
| `manual_annot_aggressive` | 0.5578 | Including drastic edits; over-corrects |
| `prot_rescue4` | 0.5578 | Only 4 Protenix swaps; undershoots |
| `combo_prot8_esm8` | 0.5563 | Protenix + ESMFold2 combined on 16 ligands |
| `manual_annot` | 0.5558 | 3-variant medchem, first run |
| `af3_iface_rescue8` | 0.5557 | AF3 interface + Protenix tail |
| `threevar_all_light` | 0.5570 | All light medchem edits |
| `rmsdpred_rescue8` | 0.5574 | Learned RMSD-predictor based selection |
| `sub_af3_chai_b2` | 0.5472 | 4-model base before Protenix layer |
| `hardtail_rescue` | 0.5456 | Tail-10 Protenix swaps |
| `hardtail20_rescue` | 0.5441 | Tail-20 Protenix swaps |
| `threevar_recommended_safe` | 0.5462 | Safe medchem swaps |
| `sub_af3_b2` | 0.5405 | 2-model |
| `prot5` | 0.5241 | First Protenix test (5 random swaps) |
| `zhyb5_prot` | 0.5241 | 5-model z-hybrid (extra model hurts) |
| `esmfold2_full` | 0.5264 | ESMFold2 alone |
| `threevar_recommended_full` | 0.5363 | Drastic medchem included |
| `combo_prot8_boltz1_8` | 0.5362 | Boltz-1 on tail |
| `sub_af3_of3` | 0.5227 | 2-model |
| `sub_af3_chai` | 0.5212 | 2-model |
| `agentic_pool_pick` | 0.5140 | 184 PXR-expert agents with SAR context — agentic reasoning cannot beat numeric selection |
| `posescorer_pick` | 0.4762 | XGBoost LambdaMART custom pose scorer — fine-ranking within good poses doesn't generalize |
| `262_ai_ipde_select` | 0.4996 | Original IPDE baseline (week 1) |

### Approaches refuted via GT-gate (never submitted)

| Approach | Refutation |
|---|---|
| OpenMM local MD refinement | 2 Å translation not recovered; closes wrong minima |
| Rotational search at fixed centroid | Too slow for 184; doesn't reliably recover native |
| Medchem agent free ligand re-drawing | 4NY9: 3.88 → 24.63 Å; agents hallucinate unphysical poses |
| Genetic anchor+tail crossover | Oracle improves but no selector can identify the improved hybrid |
| Crystal anchor priors | Fails on PanDDA fragments (bind novel sub-pockets) |
| Pocket-engagement consensus | Gameable; geometric signals don't correlate with LDDT-PLI |
| MMFF strain gating | blend_top3 = 2.251 Å vs IPDE 2.230 Å; +0.020 worse, rejected |
| Ligand-only MMFF relax of GT | Monotonically hurts (bound conformer is legitimately strained) |

---

## Why the leader scores 0.5725

The leader (`dnan-ipd`, Apheris AI) runs the **Federated OpenFold3 Initiative** — AbbVie, J&J, BMS, and Takeda co-fine-tuned OpenFold3 on thousands of experimental protein-ligand crystal structures under a privacy-preserving federated protocol. Their model has been **trained to predict PXR-family binding modes** from proprietary crystallographic data we cannot access.

The gap is **generation quality** (the model produces better poses from the start), not selection. Our pool oracle is ~1.08 Å median, which caps selection at ~0.57–0.60 LDDT-PLI. Their oracle is probably 0.70–0.80 Å.

**Our nearest public-data analog** (post-challenge): fine-tune Protenix on the NR family from the SAIR dataset (SandboxAQ, 5.2M Boltz-1x cofolds, CC BY 4.0) — scripts ready in `scripts/finetune/`.

---

## Infrastructure: automated submission ladder

The submission pipeline is fully automated to respect the ~4h rate limit and protect our best score at deadline.

### Ladder submitter (`scripts/ladder_submit.py`)

Single authority for all submissions. Runs every 4 hours via Windows Task Scheduler. On each tick:

1. Fetches leaderboard → records realized score of any in-flight candidate
2. Settles in-flight candidates before submitting next
3. Submits the next highest-priority queued candidate when cooldown elapsed
4. When queue is drained: restores `prot_rescue8.zip` if it's not currently displayed
5. After `DEADLINE_UTC = 2026-07-01 20:00 UTC`: auto-engages deadline guard (no more test candidates, only restore best)

```powershell
# Run a manual tick
.venv\Scripts\python.exe scripts\ladder_submit.py tick

# Show queue + state
.venv\Scripts\python.exe scripts\ladder_submit.py status
```

### Deadline safety net

Two layers ensure `prot_rescue8` is the final standing score:

1. **Ladder cron** (every 4h, cron `9d74752e`): auto-engages deadline guard at 20:00 UTC 07-01
2. **OS Scheduled Task** (`pxr_auto_submit_structure`): fires at 22:30 UTC 07-01, runs `deadline_restore.py`, re-submits `prot_rescue8.zip` if not already displayed

```powershell
# Verify the OS task is armed
schtasks /Query /TN "pxr_auto_submit_structure" /V /FO LIST
```

### Submission hygiene (applied before every submit)

Every submission automatically goes through three pre-submit passes in `ladder_submit.py._submit()`:

1. `clean_scoring_fails.py` — swaps `x03063`, `x03260`, `x03264` to clean fallback poses (the only 3 server-confirmed scoring failures; strict full check over-flags false positives)
2. `add_conect.py` — injects explicit CONECT records into all 184 PDBs so server reads bond topology, not 3D geometry
3. `apply_overrides.py` — forces any hand-built override poses (e.g. idealized-geometry fallback for `x03063-1`)

---

## Reproducible rebuild of the best submission

```powershell
# Environment
$env:PYTHONIOENCODING="utf-8"; $env:PYTHONUTF8="1"

# Step 1: rebuild the 4-model z-hybrid base
.venv\Scripts\python.exe scripts\build_ai_4model_zhybrid.py
#   -> submissions/zhybrid_4model.zip (baseline 0.5472)

# Step 2: apply Protenix-v2 swaps on the 8 lowest-confidence ligands
.venv\Scripts\python.exe scripts\build_prot_rescue.py 8 prot_rescue8
#   Protenix poses come from C:/tb/protenix_best/<sid>/ (archived to O:\rclone-offload\ after run)
#   -> submissions/prot_rescue8.zip (0.5640)

# Step 3: validate
.venv\Scripts\python.exe scripts\validate_submission.py submissions\prot_rescue8.zip
#   Must report: 0 errors, 184 PDBs, all LIG resname, all graph-valid

# Step 4: submit (ladder handles this automatically, but manual override:)
.venv\Scripts\python.exe scripts\submit.py submissions\prot_rescue8.zip --confirm
```

**Note on Protenix-v2 pose availability:** The full 25-sample per-ligand Protenix run was on Explorer SLURM job 7762644. The raw poses (`C:/tb/protenix_best/`) were cleaned for disk. Archived to `O:\rclone-offload\protenix_best_184\`. To rebuild from scratch: re-run Protenix on Explorer with `protenix pred --input <yaml> --sample 25 --step 200 --cycle 10`.

---

## Validation setup

The 53-holo PXR validation set gates every new candidate before submission.

Composition:

1. `scripts/build_validation_set.py` built the original 35-holo set from cached PXR structures in `data/external/pxr_holo` and `data/external/pxr_crystals`.
2. Ligands were extracted from non-solvent `HETATM` records; buffers, waters, metals, tiny ligands, molecules below 150 Da, duplicate canonical SMILES, and multi-ligand ambiguity were filtered out.
3. `scripts/fix_validation_smiles.py` replaced RDKit PDB-derived SMILES with RCSB Chemical Component Dictionary canonical SMILES, then rewrote Boltz YAML inputs using the PXR FASTA.
4. `scripts/_prepare_nb16_validation.py` added 18 ChEMBL/RCSB PXR holos from `new_validation_pdbs.csv`, all with `in_structure_test=0`.
5. Final local validation matrix: 53 holos x 20 Boltz predictions = 1060 scored poses.

Scoring harness: align predicted protein to the crystal by C-alpha Kabsch superposition, transform the ligand into the crystal frame, compute RDKit symmetry-aware heavy-atom ligand RMSD, then report `1/(1+RMSD)` as an LDDT-PLI-like proxy.

```powershell
# Score 53 holo crystals (RMSD per pose per model)
.venv\Scripts\python.exe scripts\score_validation.py data\processed\validation_set\preds_nb16
#   -> data/processed/validation_set/scores.csv
#   -> data/processed/validation_set/scores_53holo.csv snapshot

# Aggregate: method_comparison.csv (mean RMSD, % under 2A, lddt_proxy)
.venv\Scripts\python.exe scripts\aggregate_validation.py

# Pre-flight gate vs prot_rescue8 (5-30% divergence band)
.venv\Scripts\python.exe scripts\validate_selection.py submissions\<NEW>.zip submissions\prot_rescue8.zip
```

Gate threshold (from empirical calibration): `mean RMSD <= 2.15 A AND divergence in [5, 30]%`

| Method | Mean RMSD (53 holos) | Leaderboard score |
|---|---|---|
| prot_rescue8 | ~1.90 A (est.; not all final sources exist in the holo gate) | **0.5640** |
| IPDE-best (baseline) | 1.938 A | 0.4996 |
| plddt-best | 1.968 A | ~0.50 |
| oracle lower bound | 1.403 A | ~0.70 theoretical |

---

## Compute backends used

| Backend | Used for | Cost |
|---|---|---|
| Explorer (NEU HPC, `ssh explorer`) | Protenix-v2 full 184 (SLURM 7762644), ESMFold2 deep harvest, Chai-1 | $0 (free cluster) |
| Modal A100 | AlphaFold3 inference (image cached, `scripts/modal_af3.py`) | ~$50 |
| Boltz API | Boltz-2.1 best-of-5/10 (184 full + 87 drug-like massive) | ~$80 of $100 credits |
| Google Colab Pro | Chai-1 CPU fallback, Protenix testing | ~10 CU |
| dargason/cofold HF dataset | 3,680 OF3 + 3,680 B2 poses pre-built | Free (CC BY 4.0) |

---

## Post-challenge paths (not active during competition)

### 1. Protenix LoRA fine-tune on SAIR NR family (Explorer, $0)

Fine-tune Protenix on NR-weighted structures from the SAIR dataset to replicate (in public-data form) what Apheris did with proprietary pharma crystals.

```bash
ssh explorer
sbatch scripts/finetune/protenix_lora.sbatch pr_fxr_pxr
python scripts/finetune/eval_ft_on_holos.py <ckpt>   # THE GATE
```

See `scripts/finetune/README.md` for full run order. Expected cost: 1–3 days on 1×A100.

### 2. AF3-selective via Latch Bio credits (~$500)

Run AlphaFold3 with selective/custom MSA on all 184 ligands in parallel. The MSA depth is the strongest untested lever for AF3 (ESMFold2 community finding: no-MSA dramatically improves PXR score). Latch pods provide the parallelism to sweep MSA variants that Modal's quota couldn't handle.

### 3. Re-score with post-challenge tools

Once the ground-truth poses are released by OpenADMET, use them to:
- Validate the pool oracle estimates
- Train a proper PXR pose scorer (see `scripts/pxr_pose_score.py` for the XGBoost scaffold; the problem was training label noise without GT)
- Calibrate ESMFold2 MSA-depth findings vs the actual leaderboard metric

---

## Key PXR structural knowledge (used throughout)

**Pocket residues (key anchors):**

| Residue | Role |
|---|---|
| Ser247 (internal: 106) | Polar anchor, H-bond donor/acceptor |
| Gln285 (internal: 144) | Polar anchor |
| His407 (internal: 266) | Polar anchor, can be protonated |
| Arg410 (internal: 269) | Salt bridge / polar anchor |
| Phe288, Trp299, Tyr306 | Hydrophobic/aromatic subpocket |

**Drug-like vs fragment split (the key score insight):**

The public leaderboard during the challenge reflected fragment compounds only. The true final metric covers ALL 184 ligands. Drug-like compounds (MW≥330 or rotB≥5 or HA≥24; ~87/184) score ~0.46 per model — **the drug-like half is where points are left on the table**, not the fragments where all models perform similarly well at ~0.55–0.57.

**Selection wall:**

With a pool oracle of ~1.08 Å and no ground truth to train a selector, "chemically plausible ≠ native" — every smart selector we tried (learned scorer, agentic, edge-chemistry, consensus) regressed vs simple z-hybrid. The only effective lever: (1) raise the generation floor (more samples, better model), (2) trust the confidence signals as a relative ranking within each model.

---

## Environment

```
Python: .venv/Scripts/python.exe (Windows, uv-managed)
Set: PYTHONIOENCODING=utf-8 PYTHONUTF8=1
Key deps: numpy, rdkit, mdanalysis, gradio_client, scikit-learn, xgboost, scipy, openmm

Credentials (.env, gitignored):
  BOLTZ_API_KEY=...
  PXR_HF_USERNAME=xX-its-amit-Xx
  PXR_USER_ALIAS=...
  PXR_PARTICIPANT_NAME=...
  PXR_DISCORD=...
  PXR_EMAIL=...
  PXR_AFFILIATION=...
  PXR_MODEL_TAG=...
```

---

## Documentation index

- **[docs/INDEX.md](docs/INDEX.md)** — master pointer to all figures, stats, memos
- **[docs/final_submission_pixel_report.html](docs/final_submission_pixel_report.html)** — polished final submission page in pixel/late-night style
- **[docs/methodology.md](docs/methodology.md)** — IPDE selection method + validation pipeline + gates (pre-Protenix era)
- **[docs/iterations.md](docs/iterations.md)** — Iter 0–7 narrative
- **[docs/stats/methods_table.md](docs/stats/methods_table.md)** — per-method audit
- **[docs/stats/leaderboard_table.md](docs/stats/leaderboard_table.md)** — field shift table
- **[docs/stats/iter_summary.md](docs/stats/iter_summary.md)** — per-iter accounting
- **[scripts/finetune/README.md](scripts/finetune/README.md)** — post-challenge fine-tune run order
- **[submissions/af3_batches/UPLOAD_INSTRUCTIONS.md](submissions/af3_batches/UPLOAD_INSTRUCTIONS.md)** — manual AF3 server upload flow
- **[data/processed/manager_log.md](data/processed/manager_log.md)** — running action log for this session
