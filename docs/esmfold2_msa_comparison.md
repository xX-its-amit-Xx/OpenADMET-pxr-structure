# ESMFold2 / MSA-depth sweep — PXR structure challenge (leaderboard side-by-side)

**Updated:** 2026-06-24 · **Metric:** HuggingFace leaderboard LDDT-PLI (primary), with BiSyRMSD (lower=better)
and LDDT-LP (pocket LDDT) as secondary columns. **Our best overall:** `prot_rescue8` = 0.5640 (rank 2; a
4-model ensemble, not ESMFold). **Leader:** ver228 = 0.5725.

## What this experiment is
The Biohub **ESMFold2** cofolder (SMILES + protein, native MSA-free LM backbone) run as an **MSA-depth sweep**
on the 184 test ligands, to test the community/our-own hypothesis: *with less MSA the ligand drives the pose*
(less protein-template bias) → potentially better poses for the **novel PanDDA fragments** that dominate the
hard tail. Operationalized as full-MSA / shallow-MSA / no-MSA cofold variants, each a full 184-pose submission.

## Results (leaderboard)
| variant | MSA depth | LDDT-PLI | BiSyRMSD (Å) | LDDT-LP | rank | status |
|---|---|---|---|---|---|---|
| **esmfold2_full** | full MSA | **0.5264** | 3.838 | 0.9087 | 7 | ✅ scored |
| esmfold2_none | no MSA (single-seq) | _pending_ | — | — | — | ⏳ queued (prio 11.587) |
| esmfold2_shallow | shallow MSA (~8–32 seq) | _pending_ | — | — | — | ⏳ queued (prio 11.588) |
| esmfold2_rescue8 | full MSA, on 8 failure-tail only | 0.5610 | 3.483 | 0.9116 | 2 | ✅ scored (ensemble swap, not standalone) |

**Reference points:** 4-model base ensemble = 0.5551 · best ensemble `prot_rescue8` = 0.5640 · leader = 0.5725.

## Read so far (will firm up once no-/shallow-MSA score)
- **Full-MSA ESMFold2 as a standalone model is weak (0.5264)** — comparable to other single-model standalones
  (~0.52), and ~0.03 below the 4-model ensemble. As a *standalone* it does not compete; its value (if any) is
  as an ensemble member / failure-tail rescuer.
- **As a failure-tail rescuer (`esmfold2_rescue8` = 0.561)** it lands below our best (0.5640) — i.e. ESMFold2
  poses do **not** beat Protenix on the hard tail. So with **full MSA**, ESMFold2 is neutral-to-negative for us.
- The **open question this experiment answers** is whether **dropping/shrinking the MSA** flips that — if the
  no-MSA / shallow-MSA variants beat the full-MSA 0.5264 (and ideally the 0.5551 base), the "ligand-drives-pose"
  hypothesis holds and shallow-MSA becomes a real lever. **Those two numbers are pending** (jobs queued; ETA ~1–2
  days on the 1-submission-per-4h ladder). This doc auto-updates when they land.

## Caveats for the team
- **LDDT-LP is decoupled from LDDT-PLI** in our data (Spearman ≈ +0.01 across 18 submissions) — do **not** rank
  ESMFold variants by LDDT-LP; it tracks pocket-backbone quality, not ligand placement. **BiSyRMSD ≈ LDDT-PLI**
  (ρ ≈ +0.94), so it's a safe corroborating column.
- `esmfold2_rescue8` is an *ensemble swap* (ESMFold poses on 8 ligands over a 4-model base), not a pure ESMFold
  submission — included for context, not as an apples-to-apples MSA-depth cell.
- A parallel **Boltz-2 no-MSA** 184-cofold (different model, same hypothesis) was launched 2026-06-24 via the
  Boltz hosted API — a second, independent read on the MSA-depth lever; results will be added as a companion row.

_Source data: `data/processed/submission_metrics.csv` (consolidated per-submission metrics) + `data/processed/ladder_queue.csv`._
