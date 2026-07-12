# Documentation Index

Master pointer for the OpenADMET PXR Structure-track close-out. All paths relative to repo root `d:\Users\ashenoy00000\.windsurf\OpenADMET-pxr-structure\`.

## Top-level

- **[../README.md](../README.md)** — entry point: status, repo layout, current best, how to run, ceiling.
- **[final_submission_pixel_report.html](final_submission_pixel_report.html)** — polished single-page final submission narrative in pixel/late-night style.
- **[posthoc_analysis.html](posthoc_analysis.html)** + **[compound_gallery.html](compound_gallery.html)** (all-184 gallery) + **[posthoc_animation.html](posthoc_animation.html)** (3D, multi-holo) — 🔬 **post-hoc autopsy** (continually updated): per-compound traceable diagnosis for all 184 with an interactive 3D cross-model pose scatter, the quantified selection wall (confidence≠consensus in 171/184), the crystal-frame correction (boltz1/decaf frame-incompatible), three GT holo deep-dives, and rebuilt methods. Research provenance: [why_we_missed_research.md](why_we_missed_research.md).
- **[../linkedin.md](../linkedin.md)** — ready-to-post challenge wrap-up (organizers + LatchBio thanks, creative approaches, surprising findings, live links).
- **[methodology.md](methodology.md)** — full write-up of the `#262 IPDE-select` method, validation pipeline, pre-flight gates.
- **[iterations.md](iterations.md)** — iter 1–7 narrative: what was tried, what worked, what failed, why.

## Stats tables

| Table | Path | Contents |
|---|---|---|
| Methods table | [stats/methods_table.md](stats/methods_table.md) | 47 methods across iter 0–7 with verdict + reason |
| Leaderboard table | [stats/leaderboard_table.md](stats/leaderboard_table.md) | Top-20 field shift 2026-06-02 → 2026-06-08, self-overwrite casualties |
| Iter summary | [stats/iter_summary.md](stats/iter_summary.md) | Per-iter accounting; final honest assessment + ceiling |

## Figures

| # | File | Caption |
|---|---|---|
| 1 | [figures/fig01_leaderboard_trajectory.png](figures/fig01_leaderboard_trajectory.png) | Our trajectory + top-5 competitor traces (2026-05-18 → 2026-06-08) |
| 2 | [figures/fig02_validation_methods.png](figures/fig02_validation_methods.png) | Mean per-pose RMSD bar chart, 8 selectors on n=35 PXR holos |
| 3 | [figures/fig03_predicted_vs_actual.png](figures/fig03_predicted_vs_actual.png) | Predicted-vs-actual LDDT-PLI scatter (n=3 candidates with prior estimate) |
| 4 | [figures/fig04_iter_attempts.png](figures/fig04_iter_attempts.png) | 31 attempted approaches grouped by iter (status + lift) |
| 5 | [figures/iter4_validation_expansion_35_vs_53.png](figures/iter4_validation_expansion_35_vs_53.png) | Method comparison on 35 vs 53 holos (plddt overfit unmasked) |
| 6 | [figures/iter4_oracle_rmsd_per_holo.png](figures/iter4_oracle_rmsd_per_holo.png) | Per-holo oracle RMSD (lower envelope across 700+ poses) |
| 7 | [figures/iter4_kernel_closeout.png](figures/iter4_kernel_closeout.png) | nb10/15/16/17 final states (3 crashes, 1 success) |

## Memory memos (under `D:\Users\ashenoy00000\.claude\projects\d--Users-ashenoy00000--windsurf-OpenADMET-pxr-structure\memory\`)

Project-context memos (auto-persisted across sessions). Indexed by `MEMORY.md`.

### Challenge / state
- `project_pxr_structure.md` — challenge facts: 184 ligands, LDDT-PLI, P2₁2₁2₁ apo PanDDA, GitHub re-refined PDBs.
- `project_leaderboard_state.md` — 2026-06-02 snapshot; superseded by iter4_final.
- `project_omnibus_iter4_final.md` — iter4 close-out: 3 kernel crashes, 53-holo expansion, IPDE robust.
- `project_numerical_reassessment_2026_06_08.md` — gap to top-5 +0.0119, ceiling 0.515–0.535.

### Method results
- `project_validation_pipeline.md` — 35-holo validation pipeline.
- `project_validation_results_v1.md` — methods within 0.05 Å on n=35 (within noise).
- `project_validation_results_v2.md` — 5 untested methods, none beat plddt.
- `project_validation_results_v3.md` — pocket_plddt strongest at 2.0995 (-0.09 Å vs plddt).
- `project_validation_results_v3_realized.md` — hybrid_pocket_plddt scored 0.4702 (sign-inverted miss).
- `project_approach_comparison.md` — pure IPDE beats all blended approaches.
- `project_signal_correlations_v1.md` — kernel-internal signals are noise (|ρ|<0.08).
- `project_mcs_series_explored.md` — MCS series 178→120 mostly artifact; +0.008 sub-noise.
- `project_activity_structure_coupling.md` — drug-like analogs are PXR-INACTIVE; whole 184 is weak-binder.
- `project_geometric_signals_refuted.md` — GT harness; series/pocket/consensus all fail.

### Closures (this iter)
- `project_omnibus_iter1_closure.md` — 5 external paths CLOSED_NEG / BLOCKED_INFRA.
- `project_omnibus_iter3_closure.md` — iter3 pivot to external signals.
- `project_omnibus_iter4_status.md` → `project_omnibus_iter4_final.md` — kernel close-out.
- `project_mmff_strain_closure.md` — #302 MMFF strain validator FAIL; QTLSP transitively REJECTED.
- `project_ampt_closure.md` — AMPT orthogonality gate barely passed (|ρ|=0.294); GO_BUILD eligible but HOLD.
- `project_omnibus_iter5_brainstorm.md` — QTLSP GO_BUILD (now rejected by MMFF closure).
- `project_omnibus_iter6_brainstorm.md` — AMPT GO_RESEARCH, 4 ideas truncated.
- `project_documentation_v1.md` — this iter (iter 7) closeout (you are here).

### Feedback / infra
- `feedback_submission_format.md` — zip 184 PDBs, LIG resname, ligand-graph match.
- `feedback_leaderboard_shows_latest.md` — every submission overwrites; no-restore-reflex.
- `feedback_check_leaderboard_before_submit.md` — log local vs leaderboard divergence.
- `feedback_gradio_submit_api.md` — 12-arg /submit_predictions signature.
- `feedback_kaggle_p100_cuda.md` — P100 cc<7.0; T4/L4/A100 only for folding models.
- `feedback_kaggle_push_path.md` — full venv path for kaggle.exe.
- `feedback_kaggle_notebook_issues.md` — non-ASCII, rdkit-pypi dead, PDBQT spec.
- `feedback_resource_limits.md` — psutil before parallel/heavy.
- `feedback_memory_constraints.md` — RAM ceiling (93% has been hit).
- `feedback_autonomous_decisions.md` — no AskUserQuestion in /loop.

### Reference data
- `reference_dargason_cofold.md` — HF dataset: 3,680 Boltz-2 + 3,680 OF3 PDBs.
- `reference_pxr_holo_pdbs.md` — 9 PXR LBDs cached at data/external/pxr_holo/.
- `reference_pxr_pocket_map.md` — 28-residue pocket map (PXR LBD).
- `reference_external_datasets.md` — PoseBusters, Plinder, PoseBench; failed search patterns.
- `reference_sibling_repo.md` — sibling repo whitelist; reuse via sys.path.

## External resources (manual / pending)

- **AF3 batches** — [../submissions/af3_batches/UPLOAD_INSTRUCTIONS.md](../submissions/af3_batches/UPLOAD_INSTRUCTIONS.md) (10 JSONs, 184 jobs, ~10 days).
- **Re-refined PDBs** — https://github.com/OpenADMET/pxr_xtal_re-refinement (64 PDBs; HF NOT updated).
- **Re-refinement blog** — https://openadmet.ghost.io/pregnane-x-receptor-pdb-structure-rerefinement/

## Quick-start commands

```powershell
# Rebuild #262 from cache
.venv\Scripts\python.exe scripts\build_ai_ipde_select.py

# Validate
.venv\Scripts\python.exe scripts\validate_submission.py submissions\262_ai_ipde_select.zip

# Score on 53-holo validation set
.venv\Scripts\python.exe scripts\score_validation.py
.venv\Scripts\python.exe scripts\aggregate_validation.py

# Submit (rate-limited 1 / 4 h)
.venv\Scripts\python.exe scripts\submit.py submissions\262_ai_ipde_select.zip --confirm
```
