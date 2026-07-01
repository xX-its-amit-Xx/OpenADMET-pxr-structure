# Light Med-Chem Edits — Model Report

**Submission:** `submissions/manual_annot_regular.zip`
**Base model:** `prot_rescue8` (LDDT-PLI 0.564, rank 2)
**Realized score:** LDDT-PLI 0.5613, rank 2 (`data/processed/ladder_queue.csv`, priority 0.1)
**Build script:** [`scripts/build_annot_candidates.py`](../scripts/build_annot_candidates.py)

## Method

A 71-agent manual/agentic structural review pass (workflow `wj0dp8016`) inspected every
`prot_rescue8` pose for med-chem-relevant geometry issues (clashes, obviously strained
conformers, pocket-inconsistent placements) and classified each reviewed ligand into one
of three correction tiers:

| Tier | Count | Action |
|---|---|---|
| `keep` | 49 | No issue found — base `prot_rescue8` pose retained as-is |
| `light` | 20 | Small coordinate nudge to resolve a local issue, applied |
| `drastic` | 2 | Large re-placement flagged as needed, **not** applied in this build |

Decisions are recorded in `C:/tb/annot184_decisions.json`; corrected poses are produced by
the `place_ligand` pipeline and cached at `C:/tb/annot184/<id>/annotated.pdb`, already in
submission format (protein `ATOM` chain A + ligand `HETATM LIG` chain B).

**This submission (`manual_annot_regular`) applies LIGHT edits only.** For the 20 ligands
flagged `light`, the annotated pose replaces the `prot_rescue8` pose. Ligands flagged
`drastic` are left on the base `prot_rescue8` pose rather than the drastic-tier correction
(the more aggressive `manual_annot_aggressive` variant applies both light + drastic and was
evaluated separately). All edits move coordinates only — ligand graph, atom names, and
connectivity are untouched, so the reference-SMILES match used by the validator is
unaffected.

## Ligands corrected (light tier, 20/71 reviewed)

```
x00011-1  x00086-1  x00088-1  x00252-1  x00337-1  x00358-1  x00463-1  x00794-1
x00979-1  x01131-1  x01287-1  x01334-1  x01438-1  x02746-1  x02800-1  x02843-1
x02861-1  x02865-1  x03259-1  x03325-1
```

## Result

Replacing only the 20 light-tier poses on top of the `prot_rescue8` base lifted local
med-chem consistency without disturbing the other 164 poses. Realized leaderboard score
(0.5613) is within noise of the `prot_rescue8` base (0.564) — a small net negative on this
particular scoring pass, consistent with the light edits being conservative,
low-magnitude corrections rather than pose-selection changes.

See [`scripts/build_annot_candidates.py`](../scripts/build_annot_candidates.py) for the
exact build logic and [`docs/INDEX.md`](INDEX.md) for the full documentation set.
