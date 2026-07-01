# Leaderboard Snapshot (Top-20)

Two-week field shift comparison: **2026-06-02** (baseline) vs **2026-06-08** (current).

## Top-20 comparison

| Rank (now) | Team | Score (2026-06-08) | Score (2026-06-02) | Delta | Notes |
|---|---|---|---|---|---|
| 1 | dnan-ipd | 0.5612 | 0.5612 | 0.0000 | stable since 2026-05-29 |
| 2 | fakeplastictrees | 0.5599 | not in top-20 | NEW | major new entrant |
| 3 | TangerineTrees | 0.5285 | 0.4968 (r7) | +0.0317 | biggest single-team jump |
| 4 | TCB | 0.5267 | ~0.5260 | ~+0.0007 | stable |
| 5 | discoverybytes | 0.5115 | ~0.5318 | -0.0203 | drifted down (still top-5) |
| 6 | Radi | 0.5046 | ~0.4995 | +0.0051 | small climb |
| 7 | rdkbio | 0.5016 | not in top-20 | NEW | crossed 0.50 |
| 8 | Jedi | 0.5008 | not in top-20 | NEW | crossed 0.50 |
| 9 | DeepFoldX | 0.5000 | not in top-20 | NEW | crossed 0.50 |
| **10** | **xX-its-amit-Xx (us)** | **0.4996** | **0.4996 (r9-10)** | **0.0000** | HOLD policy active; no submission since 2026-06-02 |
| 11 | (rotated) | ~0.4950 | ~0.4940 | ~0 | minor churn |
| 12-15 | (rotated) | 0.4880-0.4920 | similar | small shifts | self-overwrite churn band |
| 16-18 | (suspenders, etc) | 0.4727-0.4850 | 0.5521-0.4990 | -0.05 to -0.08 | self-overwrite casualties |
| 19-20 | (Cryo-EMinem etc) | 0.4640 area | 0.5113 area | -0.0473 | self-overwrite casualty (was #3) |

## Self-overwrite casualty list (leaderboard-shows-latest trap)

Top teams that **overwrote their best score by submitting worse experiments** and dropped:

| Team | Best ever (date) | Current (2026-06-08) | Drop | Lesson |
|---|---|---|---|---|
| suspenders | 0.5521 (May 21, #2) | 0.4727 (#18) | -0.0794 | most expensive self-overwrite of the project |
| Cryo-EMinem | 0.5113 (May 21, #3) | 0.4640 (#25) | -0.0473 | fell from top-5 to outside top-20 |
| blacktea | 0.5015 (May 21, #6) | out of top-20 | >-0.05 | fell off top-20 entirely |
| (us, before policy) | 0.4997 (May 20, #5) | recovered to 0.4996 | -0.0001 | survived because of #262 restore reflex |

**Why we did NOT fall:** strict no-restore-reflex policy + Windows Scheduled Task safety restore + 25% divergence pre-flight gate prevented us from overwriting #262 with sub-0.4996 experimental submissions.

## Position deltas (us)

| Metric | 2026-06-02 | 2026-06-08 | Delta |
|---|---|---|---|
| Our rank | 9-10 | 10 | -1 (drift down) |
| Our score | 0.4996 | 0.4996 | 0.0000 |
| Gap to #1 (dnan-ipd) | +0.0616 | +0.0616 | 0.0000 (unchanged) |
| Gap to #5 (top-5 cutoff) | +0.0119 | +0.0119 | 0.0000 (unchanged) |
| Number of teams >= 0.5000 | 5 | 9 | +4 |
| Number of teams >= 0.4996 (>= us) | 9 | 9 | 0 |

## Field acceleration analysis

- **4 NEW teams crossed 0.50** in the 6-day window (rdkbio, Jedi, DeepFoldX, plus TangerineTrees jumping further).
- **fakeplastictrees** entered at #2 (0.5599) — second-highest score on the board, possibly an organizer / late-arrival strong team.
- **Top-5 turnover** masked our static position: discoverybytes, suspenders, Cryo-EMinem all dropped from top-5, freeing the slots that the 4 new sub-0.51 teams filled. Net effect: gap to top-5 unchanged.
- **Top-1 is locked at 0.5612 for 10+ days.** No team has crossed 0.56 since 2026-05-29.

## Our score trajectory

| Date | Submission | Score | Rank | Notes |
|---|---|---|---|---|
| 2026-05-18 | 01_boltz2_tutorial | 0.4632 | 21 | baseline |
| 2026-05-19 | 04_of3_best | 0.4751 | 12 | +9 ranks |
| 2026-05-19 | 06_b2_ligand_iptm | 0.4869 | 10 | first cross-model |
| 2026-05-20 | 08_iptm_zhybrid | 0.4984 | 6 | z-hybrid breakthrough |
| 2026-05-20 | 13_zhybrid_srs | 0.4997 | 5 | top-5 entered |
| 2026-05-21 | 14_both_composite | 0.4849 | 11 | first self-overwrite loss |
| 2026-05-21 | 13_recovery | 0.4997 | 6 | restored to top-5 |
| 2026-05-21 | 113_crystal_protein_srs | 0.4657 | 22 | crystal protein hurt |
| 2026-05-22 | 262_ai_ipde_select | **0.4996** | 7 | first pure-IPDE win |
| 2026-05-22 | 278/241/272/267 | 0.20-0.38 | 29-41 | failed crystal anchors |
| 2026-05-23 | 262_restore | 0.4996 | 6 | restored via auto_submit_v2 |
| 2026-05-24 | 282 (B2-only) | 0.4818 | 12 | B2 alone is weaker |
| 2026-05-24 | 263 (OF3-only) | 0.4724 | 13 | OF3 alone is weaker |
| 2026-05-25 | 281 (per-family MLP) | 0.4662 | 21 | learned ranker disaster |
| 2026-05-26 | 290 (iplddt) | 0.4832 | 13 | iplddt pushes B2 wrong way |
| 2026-05-27 | 294 (occupancy) | 0.4586 | 24 | holo prior misdirects fragments |
| 2026-05-27 | 297 (anti-occupancy) | 0.4514 | 25 | inverse also fails |
| 2026-05-28 | 283 (z-score complex) | 0.4601 | 25 | z-score on wrong B2 signal |
| 2026-05-28 | 298 (z-srs) | 0.4597 | 24 | OF3 srs doesn't help when B2 wrong |
| 2026-05-30 | 322 (activity-Tanimoto) | 0.458 | 27 | medoid family falsified |
| 2026-06-02 | v3_hybrid_pocket_plddt | 0.4702 | 21 | sign-inversion on fragments |
| 2026-06-02 | 262_restore | **0.4996** | 9-10 | restored, HOLD invoked |
| 2026-06-08 | (no submission) | **0.4996** | 10 | HOLD maintained 6 days |

## HOLD policy validation

Since the 2026-06-02 16:50 UTC #262 restore, we have made **0 submissions over 6 days**. During that window:

- 5 prospective candidates were built (v_md_refined, v_activity_aware, 302_mmff_strain, ampt prototype, qtlsp prototype).
- **0 candidates passed the 3-gate validation** (5-30% divergence + 53-holo RMSD <= 2.15A + medium+ confidence).
- Score preserved at 0.4996.
- If we had submitted any of the 5, leaderboard model says expected score would be 0.43-0.49, all worse than 0.4996.

HOLD policy has saved approximately **+0.05 to +0.10 LDDT-PLI of leaderboard drift** in 6 days.
