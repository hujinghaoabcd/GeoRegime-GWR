# Synthetic Simulation Status

Last updated: 2026-09-01

## Purpose

Synthetic experiments provide known ground truth for GR-GWR boundary recovery and coefficient recovery.  They complement the Georgia real-data diagnostics, where the true regimes and true local coefficients are unknown.

All results below are single-realization behavioral prototypes.  They are not yet Monte Carlo evidence and must not be reported as final performance estimates.

## Simulation 1 — canonical piecewise-constant regime discontinuity

Design:

- 25 x 25 regular lattice, n=625;
- Queen topology;
- K_true=3 contiguous regimes;
- K=3 is supplied to the estimated method, but true labels/boundary locations are hidden;
- three correlated covariates, rho=0.2;
- piecewise-constant true slopes with sharp cross-regime discontinuities;
- noise sigma=0.75;
- current estimated chain: pilot GWR -> Queen-constrained Ward -> lambda=0.5 refinement -> unified RegimeAwareGWR;
- Oracle GR-GWR uses the true labels.

Key boundary recovery:

- initial Ward ARI = 0.7671, boundary F1 = 0.2490;
- refined GR-GWR ARI = 0.9594, boundary F1 = 0.8522;
- refined precision = 0.8376, recall = 0.8673;
- true boundary edges = 113, refined estimated boundary edges = 117.

Model fit RMSE:

- OLS = 0.9535;
- GWR = 0.4552;
- MGWR = 0.4263;
- estimated GR-GWR = 0.4077;
- Oracle GR-GWR = 0.4055.

True boundary-jump recovery ratios:

- X1: GWR 0.144, MGWR 0.165, GR-GWR 0.760, Oracle 0.953;
- X2: GWR 0.178, MGWR 0.236, GR-GWR 0.852, Oracle 1.060;
- X3: GWR 0.173, MGWR 0.255, GR-GWR 0.971, Oracle 0.984.

Interpretation: the prototype behaves as intended under a clean discontinuous-regime truth.  GWR/MGWR strongly attenuate cross-boundary jumps; estimated GR-GWR recovers most of them and approaches the oracle fit.

## Simulation 2 — smooth within regimes, discontinuous between regimes

Design:

- same 25 x 25 lattice, Queen graph and K_true=3;
- true coefficients contain continuous spatial fields inside regimes plus regime-specific offsets across boundaries;
- therefore the truth is not piecewise constant;
- rho=0.2, noise sigma=0.75;
- same current estimated GR-GWR chain and Oracle comparator.

The standardized true within-regime coefficient SD is non-zero for every slope and every regime, confirming that this is genuinely a continuous-within / discontinuous-between truth.

Boundary recovery:

- initial Ward ARI = 0.8819, boundary F1 = 0.4229;
- refined GR-GWR ARI = 0.9328, boundary F1 = 0.6724;
- refined precision = 0.6555, recall = 0.6903;
- true boundary edges = 113, refined estimated boundary edges = 119;
- refinement accepts 2 sweeps and then the next sweep is rejected by the global objective guard.

Model fit RMSE:

- OLS = 0.9625;
- GWR = 0.4906;
- MGWR = 0.4856;
- estimated GR-GWR = 0.4676;
- Oracle GR-GWR = 0.4678.

Important: the estimated GR-GWR slightly beating Oracle in this single in-sample realization is not evidence of superiority; estimated boundaries are response-adaptive and can exploit sample noise.  Only out-of-sample / Monte Carlo comparisons can interpret that difference.

Coefficient RMSE:

- X1: GWR 0.2147, MGWR 0.2049, GR-GWR 0.1554, Oracle 0.0825;
- X2: GWR 0.1733, MGWR 0.1685, GR-GWR 0.1133, Oracle 0.0561;
- X3: GWR 0.1658, MGWR 0.1714, GR-GWR 0.1546, Oracle 0.0523.

True boundary-jump recovery ratios:

- X1: GWR 0.186, MGWR 0.207, GR-GWR 0.721, Oracle 1.022;
- X2: GWR 0.192, MGWR 0.205, GR-GWR 0.763, Oracle 1.055;
- X3: GWR 0.201, MGWR 0.171, GR-GWR 0.792, Oracle 1.173.

### Critical diagnostic exposed by Simulation 2

Estimated GR-GWR does not uniformly dominate GWR for recovery of very small within-regime adjacent gradients.  On true within-regime Queen edges, gradient RMSE is larger for estimated GR-GWR than for GWR/MGWR because remaining false estimated boundaries create artificial coefficient jumps inside a true regime.

This is visible in the within-regime gradient RMSE:

- X1: GWR 0.0735, MGWR 0.0770, estimated GR-GWR 0.1103, Oracle GR-GWR 0.0169;
- X2: GWR 0.0718, MGWR 0.0739, estimated GR-GWR 0.0884, Oracle GR-GWR 0.0099;
- X3: GWR 0.0665, MGWR 0.0545, estimated GR-GWR 0.1331, Oracle GR-GWR 0.0098.

The Oracle result is crucial: when the boundary is correct, regime-aware refitting does not create the same problem and is much closer to the true within-regime gradients.  The degradation is therefore primarily linked to boundary-identification error / false boundaries, not simply to the idea of regime-aware local regression itself.

This leads to an important, defensible method statement:

> GR-GWR trades cross-boundary smoothing bias against boundary-identification error.  Its advantage should increase as true discontinuities become stronger and identifiable; when boundaries are weak or ambiguous, false boundaries can introduce local coefficient discontinuities.

This statement should be tested directly in Simulation 3 (boundary strength x noise), rather than hidden by tuning Simulation 2 until every metric favors GR-GWR.

## Next synthetic experiment

Simulation 3 should vary both regime discontinuity strength and noise/SNR.  The goal is to estimate the operating region where the reduction in cross-boundary smoothing bias outweighs boundary-identification error.  This is now scientifically motivated by Simulation 2 rather than being a generic sensitivity exercise.
