# ADR-0002: The extracted GR-GWR implementation is a baseline, not the paper definition

Status: Accepted
Date: 2026-08-30

## Context

The first code placed in this repository follows the current pyGWRx GR-GWR logic so that research can start from a known working idea. However, several components are already under theoretical review: kNN+MST adjacency, coordinate-weighted clustering features, constrained Ward initialization, ICM refinement, shared-regime structure, and graph-cut boundary penalty.

## Decision

The file `src/georegime_gwr/grgwr.py` is designated `GRGWRBaseline` rather than the final `GRGWR` paper implementation.

No downstream document may treat its current algorithmic choices as frozen merely because they exist in code.

A future paper-ready model should be introduced only after targeted literature checks and ablation/simulation evidence justify the retained components.

## Consequences

- Breaking algorithmic changes are allowed.
- Baseline behavior should remain reproducible long enough to compare alternatives.
- When a candidate paper algorithm is accepted, create a new ADR defining its mathematical contract before renaming or replacing the baseline implementation.
