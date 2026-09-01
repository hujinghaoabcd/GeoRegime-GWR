# GR-GWR Algorithm Review Checkpoint

Last updated: 2026-09-01

## Current decision

Pause further work on K-selection, spatial cross-validation, 1-SE rules, and other model-selection complexity.

The next task is to review the current GR-GWR algorithm from the beginning, in plain language first and mathematical detail second, before adding any new machinery.

## Working baseline to review

1. Standard GWR pilot fit provides local coefficient vectors.
2. Geography is represented by two distinct objects:
   - metric distance matrix D for GWR kernel weighting;
   - adjacency/topology matrix W (currently Queen for polygon data) for allowed spatial merges/moves and connectivity.
3. Pilot slope coefficients are standardized and used as the feature space for initial spatially constrained Ward clustering.
4. K is currently treated as externally supplied during algorithm-behavior experiments; no final K-selection rule is frozen.
5. Initial labels are refined only through spatially adjacent candidate regimes, with source-regime connectivity protection and a minimum regime size.
6. Current exploratory refinement local cost is leave-one-out squared prediction error plus lambda times Queen-neighbor label disagreement; lambda=0.5 is the current working baseline, not a final theoretical value.
7. After accepted refinement sweeps, the model is refit using one unified RegimeAwareGWR:
   w_ij = K(d_ij / b_i) * I(z_i = z_j),
   with adaptive bandwidth geometry defined within the focal observation's regime; the current working bandwidth policy is regime_size.
8. Cross-regime weights are exactly zero in the current hard-boundary baseline.
9. K-selection, lambda selection, spatial CV, effective-complexity accounting, soft-boundary variants, and computational optimization remain open issues and are intentionally parked during the algorithm review.

## Review objective

Reconfirm what each step is doing, why it is needed, which parts are core GR-GWR ideas, and which parts are only temporary experimental choices before resuming method development.
