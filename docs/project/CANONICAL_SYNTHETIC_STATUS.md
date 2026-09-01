# Canonical synthetic experiments status

Date: 2026-09-01

## Common experimental frame

- Regular 25 x 25 lattice, n=625.
- Three contiguous true regimes; K=3 is held fixed in Simulations 1–3 so boundary detectability is isolated from K selection.
- True boundary locations/labels are hidden from estimated GR-GWR.
- X and y are globally z-standardized before fitting; true coefficients are transformed to the same fitted scale.
- Comparators: OLS where relevant, external standard GWR, external MGWR, estimated GR-GWR, oracle GR-GWR.
- Estimated GR-GWR path: pilot GWR slopes -> Queen-constrained Ward -> lambda=0.5 boundary refinement -> unified RegimeAwareGWR.
- Oracle GR-GWR uses exactly the same refit primitive with true labels.

## Simulation 1 — piecewise-constant regime truth

Initial Queen-Ward: ARI=0.7671, boundary F1=0.2490.
After six accepted refinement iterations: ARI=0.9594, precision=0.8376, recall=0.8673, boundary F1=0.8522, estimated boundary edges=117 vs 113 true.

Model RMSE: GWR=0.4552, MGWR=0.4263, GR-GWR=0.4077, oracle=0.4055.

Slope coefficient RMSE:

| term | GWR | MGWR | GR-GWR | oracle |
|---|---:|---:|---:|---:|
| X1 | 0.2282 | 0.2168 | 0.1670 | 0.0380 |
| X2 | 0.1822 | 0.1628 | 0.1371 | 0.0341 |
| X3 | 0.1642 | 0.1456 | 0.0364 | 0.0375 |

True-boundary jump recovery ratios: GWR/MGWR recover only about 14–25% of the true jump, estimated GR-GWR about 76–97%, and oracle about 95–106%.

Result directory: `results/synthetic/canonical_regime_boundary/`.

## Simulation 2 — smooth within regimes + discontinuous between regimes

This is deliberately harder and less favorable to a regime model: coefficients vary continuously within each true regime and also jump across true regime boundaries.

Boundary recovery:

- initial ARI=0.8819 -> refined ARI=0.9328;
- initial boundary F1=0.4229 -> refined F1=0.6724;
- final precision=0.6555, recall=0.6903.

Model RMSE:

- GWR=0.4906;
- MGWR=0.4856;
- estimated GR-GWR=0.4676;
- oracle GR-GWR=0.4678.

Slope coefficient RMSE:

| term | GWR | MGWR | GR-GWR | oracle |
|---|---:|---:|---:|---:|
| X1 | 0.2147 | 0.2049 | 0.1554 | 0.0825 |
| X2 | 0.1733 | 0.1685 | 0.1133 | 0.0561 |
| X3 | 0.1658 | 0.1714 | 0.1546 | 0.0523 |

True-boundary jump recovery remains much stronger for estimated GR-GWR (72–79%) than GWR/MGWR (17–21%).

Important diagnostic: estimated GR-GWR is worse than GWR/MGWR for recovery of very small *within-regime* gradients because imperfect estimated boundaries create artificial local jumps. Oracle GR-GWR recovers within-regime gradients very well. This isolates a core method trade-off: reducing cross-boundary smoothing bias versus the cost of boundary-identification error.

Result directory: `results/synthetic/smooth_within_discontinuous_between/`.

## Simulation 3 — boundary strength x noise phase diagram

Single-seed 5 x 4 prototype: boundary strength delta={0,0.5,1.0,1.5,2.0}; noise sigma={0.5,0.75,1.0,1.5}; 20 cells total. Delta=0 is retained only as a forced-K stress cell and is not an identifiable boundary-recovery target.

Across the 16 positive-boundary cells:

- estimated GR-GWR has lower slope coefficient RMSE than GWR in 15/16 cells;
- lower slope coefficient RMSE than MGWR in 15/16 cells;
- lower than the better of GWR/MGWR in 15/16 cells;
- mean refined boundary F1=0.6114;
- mean refined ARI=0.8976.

Mean behavior by boundary strength:

| delta | boundary F1 | ARI | GR-GWR minus GWR coef-RMSE | GR-GWR minus MGWR coef-RMSE | jump recovery |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.423 | 0.875 | -0.0097 | -0.0072 | 0.459 |
| 1.0 | 0.533 | 0.888 | -0.0276 | -0.0231 | 0.572 |
| 1.5 | 0.752 | 0.932 | -0.0527 | -0.0454 | 0.810 |
| 2.0 | 0.738 | 0.895 | -0.0331 | -0.0220 | 0.815 |

The strongest average advantage occurs around delta=1.5 in this single seed. Delta=2.0 does not monotonically improve everything because boundary refinement is data-adaptive and local-optimum/order effects remain.

Mean behavior by noise:

| noise sigma | boundary F1 | ARI | GR-GWR minus GWR coef-RMSE | GR-GWR minus MGWR coef-RMSE |
|---:|---:|---:|---:|---:|
| 0.5 | 0.707 | 0.912 | -0.0352 | -0.0260 |
| 0.75 | 0.573 | 0.885 | -0.0265 | -0.0198 |
| 1.0 | 0.656 | 0.923 | -0.0446 | -0.0391 |
| 1.5 | 0.509 | 0.870 | -0.0168 | -0.0128 |

The one positive-boundary loss is the weakest-boundary/highest-noise cell (delta=0.5, sigma=1.5): GR-GWR slope coefficient RMSE=0.0999 versus GWR=0.0978 and MGWR=0.0954. This is scientifically useful: when discontinuity is weak relative to noise, hard boundary errors can outweigh the benefit of preventing cross-boundary smoothing.

At delta=0, forcing K=3 is clearly harmful for coefficient recovery: GR-GWR is worse than GWR/MGWR in all four no-discontinuity stress cells. This is expected and reinforces that formal K selection / K=1 negative controls are essential; it must not be interpreted as a valid test of the final method because K is intentionally fixed incorrectly.

Result directory: `results/synthetic/boundary_noise_phase/`.

## Current interpretation

The synthetic evidence now supports a coherent mechanism rather than a single favorable example:

1. GWR/MGWR can strongly attenuate known discontinuous coefficient jumps.
2. Estimated GR-GWR can recover much of those jumps without being given the true boundary locations.
3. It still works when coefficients vary smoothly inside regimes, but boundary-label errors can create artificial local discontinuities.
4. Benefits strengthen as discontinuity becomes identifiable relative to noise; weak-boundary/high-noise cases are the expected failure zone.
5. Hard masking is therefore only defensible together with reliable regime identification and a mechanism that can avoid forcing regimes when none are supported.

## Important limitations

- Simulations 1–3 are still single-realization / single-seed prototypes, not publication-grade Monte Carlo evidence.
- K=3 is supplied in all three; K selection remains OPEN.
- lambda=0.5 is a working value, not theoretically selected.
- Current refinement is sequential and can reach local optima.
- Current refit bandwidth policy is still the working `regime_size` baseline.
- No spatial cross-validation or independent test realization is used here; coefficient truth recovery is the primary validation target.
- No irregular real-world topology experiment yet.

## Next simulation tasks

1. Simulation 4: Georgia-topology semi-synthetic experiment with known irregular connected regimes.
2. Simulation 5: formal no-boundary controls (stationary and smooth-continuous truth) coupled to candidate K including K=1; this must test false-regime creation rather than force K=3.
3. Then convert the finalized scenarios to paired Monte Carlo repetitions and confidence intervals.
4. Add variable-specific-boundary stress test and computational scaling as supplementary experiments.
