# GR-GWR Refinement Computational Complexity Notes

Last updated: 2026-08-31

## Status

This note records an OPEN implementation/method issue discovered after the current Georgia K=6 unified label-refinement experiment. It is not a frozen design decision.

## Current exploratory refinement cost

The current refinement algorithm does **not** run a complete GWR for every candidate label. For each spatial unit and each admissible candidate regime, it performs a temporary **local weighted regression** to obtain a leave-one-out prediction cost. Candidate regimes are restricted to the current regime plus regimes present among Queen-neighboring units.

After one full sweep of candidate moves, the algorithm performs one complete unified `RegimeAwareGWR` refit and applies the global objective guard.

Therefore the main computational work is approximately:

`number of refinement sweeps × number of checked units × number of admissible neighboring regimes × local weighted-regression cost`

plus one complete unified `RegimeAwareGWR` refit per accepted/proposed sweep.

For the current Georgia experiment (159 counties, K=6, 8 accepted refinement iterations), this is fully tractable. However, the same brute-force pattern may become expensive for datasets with thousands or tens of thousands of spatial units.

## Important distinction: computational vs statistical complexity

The large number of temporary local regressions in the search process does **not** mean the final model has thousands of effective degrees of freedom.

- **Computational complexity**: high because many temporary local regressions are evaluated during label search.
- **Statistical model complexity**: characterized separately (for the current fitted baseline by quantities such as `trace(S)`), and is not multiplied directly by the number of candidate evaluations.

The current K=6 unified baseline has `trace(S)` around 39.7 before refinement and around 40.7 after the exploratory refinement. The search procedure itself nevertheless adds data-adaptive partition-selection complexity that is not represented by the naive fixed-partition `trace(S)`/AICc calculation. This remains an OPEN inference/model-selection issue.

## Why the current implementation is only an exploratory baseline

The current implementation sweeps over all units even though many interior units have no alternative neighboring regime and therefore cannot change label. This is acceptable for validating algorithm behavior on Georgia but should not be treated as the final scalable implementation.

## Required future efficiency work

Potential optimizations to evaluate later:

1. **Boundary-unit screening**: only evaluate units adjacent to at least one different regime; interior units cannot move under the current topology rule.
2. **Active-set refinement**: after a move, reconsider only the moved unit and nearby units whose candidate set or local objective can change.
3. **Distance/kernel caching**: reuse pairwise distances, sorted regime-local neighbor orders, and kernel components instead of recomputing them.
4. **Sufficient-statistic / matrix caching**: reuse weighted cross-products where possible and update only affected regime-local terms.
5. **Incremental refit**: avoid rebuilding every local regression when only a small part of the partition changes, if a mathematically equivalent update can be derived.
6. **Parallel candidate evaluation**: evaluate independent boundary candidates in parallel where sequential-update semantics are not required.
7. **Benchmark scaling**: explicitly measure runtime/memory against n and number of boundary units before declaring the algorithm scalable.

## Current decision

Do **not** redesign the method around computational efficiency yet. Keep the current unified `RegimeAwareGWR` + Queen-constrained refinement as the behavioral baseline while the method is still being studied. Record computational scaling as a mandatory issue to solve before the final software/algorithm specification is frozen.
