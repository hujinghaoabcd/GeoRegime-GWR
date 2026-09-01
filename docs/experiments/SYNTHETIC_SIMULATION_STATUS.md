# Synthetic Simulation Status

Last updated: 2026-09-01

## Scope and evidence level

Synthetic experiments provide known ground truth for GR-GWR boundary and coefficient recovery. They complement the Georgia real-data diagnostics, where the true regimes and true local coefficients are unknown.

**Current evidence level:** Simulations 1–4 are behavioral prototypes. Simulation 3 covers a 5 x 4 parameter grid but still uses one seed per cell. None of the results below should yet be presented as final Monte Carlo performance estimates. K is fixed to the known experimental K in the current prototypes so that boundary recovery can be studied separately from K selection.

Current estimated chain used throughout:

`pilot GWR -> Queen-constrained Ward -> lambda=0.5 boundary refinement -> unified RegimeAwareGWR`

Oracle GR-GWR uses the same final refit primitive but receives the true regime labels.

## Simulation 1 — canonical piecewise-constant discontinuity

Design: 25 x 25 lattice (n=625), Queen topology, K_true=3, rho_X=0.2, noise sigma=0.75. True slopes are piecewise constant with sharp regime boundaries; true labels are hidden from estimated GR-GWR.

Key results:

- initial ARI 0.7671 -> refined ARI 0.9594;
- boundary F1 0.2490 -> 0.8522;
- fit RMSE: GWR 0.4552, MGWR 0.4263, GR-GWR 0.4077, Oracle 0.4055;
- jump recovery: X1 GWR 0.144 / MGWR 0.165 / GR-GWR 0.760 / Oracle 0.953;
- X2: 0.178 / 0.236 / 0.852 / 1.060;
- X3: 0.173 / 0.255 / 0.971 / 0.984.

Interpretation: under clean discontinuous truth, GWR/MGWR strongly attenuate true jumps while estimated GR-GWR recovers most of them and approaches the oracle fit.

## Simulation 2 — smooth within, discontinuous between

Design: same lattice and K_true=3, but every coefficient contains genuine smooth within-regime variation plus regime-specific discontinuous offsets.

Key results:

- initial ARI 0.8819 -> refined ARI 0.9328;
- boundary F1 0.4229 -> 0.6724;
- fit RMSE: GWR 0.4906, MGWR 0.4856, GR-GWR 0.4676, Oracle 0.4678;
- slope coefficient RMSE: X1 GWR 0.2147 / MGWR 0.2049 / GR 0.1554 / Oracle 0.0825;
- X2: 0.1733 / 0.1685 / 0.1133 / 0.0561;
- X3: 0.1658 / 0.1714 / 0.1546 / 0.0523;
- true jump recovery: X1 0.186 / 0.207 / 0.721 / 1.022; X2 0.192 / 0.205 / 0.763 / 1.055; X3 0.201 / 0.171 / 0.792 / 1.173.

Critical diagnostic: on true within-regime Queen edges, estimated GR-GWR has worse tiny-gradient recovery than GWR/MGWR when false inferred boundaries remain. Oracle GR-GWR does not show this problem. Therefore the main trade-off is:

> GR-GWR reduces cross-boundary smoothing bias but pays a cost for boundary-identification error. False boundaries can create artificial local discontinuities.

This adverse result is retained deliberately and motivates Simulation 3.

## Simulation 3 — boundary strength x noise phase diagram

Two independently implemented single-seed 5 x 4 grids have been run. The more structured scenario-matrix version uses delta = 0, 0.5, 1.0, 1.5, 2.0 and noise sigma = 0.5, 0.75, 1.0, 1.5, with one common seed per cell and K fixed to 3.

Across the 16 cells with a genuine positive boundary contrast:

- GR-GWR slope coefficient RMSE beats GWR in 15/16 cells;
- beats MGWR in 15/16 cells;
- beats the better of GWR/MGWR in 15/16 cells;
- mean refined ARI = 0.8976;
- mean boundary F1 = 0.6114.

By boundary strength, mean GR-GWR minus GWR coefficient RMSE (negative is better):

- delta=0.5: -0.0097, mean F1=0.4228, mean jump recovery=0.4593;
- delta=1.0: -0.0276, F1=0.5328, jump=0.5719;
- delta=1.5: -0.0527, F1=0.7521, jump=0.8095;
- delta=2.0: -0.0331, F1=0.7377, jump=0.8148.

The pattern is broadly consistent with stronger true discontinuities being easier and more valuable to model, but it is not monotonic in every single cell because this is still one realization and the sequential refinement can enter different local optima.

### Critical null-boundary warning

At delta=0 there is no statistically identifiable coefficient boundary. K=3 is nevertheless forced in this prototype. GR-GWR can slightly improve **in-sample fitted RMSE** while making coefficient recovery substantially worse than GWR/MGWR. For example at noise sigma=0.75:

- GWR coefficient RMSE = 0.02156;
- MGWR = 0.02199;
- forced-K GR-GWR = 0.05040;
- yet GR-GWR in-sample fit RMSE is slightly lower.

This is direct evidence that in-sample fit cannot be used to justify segmentation. Delta=0 ARI/F1 is intentionally treated as undefined rather than a valid recovery target.

Publication-level Simulation 3 therefore requires paired Monte Carlo replications and an explicit no-boundary/K-selection control.

## Simulation 4 — Georgia real-topology semi-synthetic experiment

Design:

- real Georgia 159-county geometry and UTM coordinates;
- real verified Queen graph with 431 edges;
- K_true=4;
- true regimes generated independently of y by deterministic multi-source Queen-graph growth from four geographically separated seed counties, so every true regime is connected and boundaries are irregular;
- true regime sizes = 39, 52, 29, 39;
- synthetic X and y with smooth within-regime coefficient variation plus discrete regime shifts;
- rho_X=0.2, noise sigma=0.75;
- K=4 supplied, true labels/boundaries hidden from estimated GR-GWR.

Boundary recovery:

- initial Queen-Ward ARI = 0.5794;
- refined ARI = 0.7858;
- initial boundary F1 = 0.2481;
- refined boundary F1 = 0.6769;
- true boundary edges = 62;
- final estimated boundary edges = 68;
- precision = 0.6471, recall = 0.7097;
- 3 refinement sweeps accepted, stop reason `no_label_changes`.

Model recovery:

- GWR: fit RMSE 0.6266, slope coefficient RMSE 0.2910, true-jump recovery 0.1638;
- MGWR: 0.5362, 0.2956, 0.2054;
- estimated GR-GWR: 0.4842, 0.2394, 0.6612;
- Oracle GR-GWR: 0.4606, 0.1016, 0.9977.

Interpretation: the core behavior survives replacement of the regular lattice by real irregular county geometry, metric distances and Queen topology. Boundary recovery is harder than the clean lattice experiment, but refinement materially improves the initial spatial partition and GR-GWR recovers much more of the true coefficient jump than GWR/MGWR. The remaining gap to Oracle again shows that boundary-identification error, not only the final regime-aware refit, is the central unresolved difficulty.

## Current conclusions from Simulations 1–4

1. The intended discontinuity mechanism works under known ground truth.
2. The benefit is not limited to piecewise-constant coefficient surfaces.
3. Stronger identifiable boundaries generally increase the value of GR-GWR.
4. False boundaries are costly; forced segmentation can overfit even when there is no true boundary.
5. The behavior persists on real irregular Georgia geography/topology.
6. Current single-run evidence is not enough for final superiority claims.

## Next synthetic work

Highest priority before publication-level claims:

1. paired multi-seed / Monte Carlo replication of representative Simulation 3 cells;
2. explicit no-boundary negative control without forcing K>1;
3. K-selection behavior;
4. then variable-specific boundary mismatch and computational scaling sensitivity.
