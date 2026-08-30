# ADR-0001: Research repository scope

Status: Accepted
Date: 2026-08-30

## Context

GR-GWR originally lived inside pyGWRx, whose architecture is designed for a broader, more mature software package. Method-paper development needs a much smaller environment where the algorithm can change rapidly, including changes that may invalidate the current class structure or optimization procedure.

## Decision

Create an independent repository, `GeoRegime-GWR`, with the following rules:

1. Keep only a minimal baseline GWR and a readable GR-GWR research implementation.
2. Treat the current GR-GWR implementation as a baseline snapshot, not a frozen API or final model.
3. Optimize for mathematical transparency, experiments, ablations, and reproducibility rather than package completeness.
4. Allow major algorithm components to be deleted or replaced.
5. Record every paper-defining design change in ADRs and `CURRENT_STATUS.md`.
6. Maintain `00_PROJECT_HANDOFF.md` as the cross-conversation recovery entry point.

## Consequences

Positive:

- GR-GWR can evolve independently from pyGWRx compatibility constraints.
- Algorithmic alternatives can be compared without destabilizing the mature package.
- New conversations can recover the exact research state from repository documents.

Negative:

- Some code is duplicated from the conceptual baseline in pyGWRx.
- The research repository is not intended to provide production-grade API stability.
- Improvements made here do not automatically propagate back to pyGWRx.

## Follow-up

Only after the paper algorithm is stabilized should a deliberate decision be made about whether and how to port the final method back into pyGWRx.
