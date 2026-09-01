# Canonical synthetic prototype status

Date: 2026-09-01

## Purpose

Single-run behavioral validation for GR-GWR under a known piecewise-constant spatial coefficient truth.

- 25 x 25 regular lattice, n=625.
- Three contiguous true regimes.
- K=3 is fixed for this first prototype, but true boundary locations and labels are hidden from estimated GR-GWR.
- X and y are globally z-standardized before fitting.
- Compared models: OLS, standard GWR, MGWR, estimated GR-GWR, oracle GR-GWR.
- Estimated GR-GWR path: pilot GWR local slopes -> Queen-constrained Ward K=3 -> lambda=0.5 boundary refinement -> unified RegimeAwareGWR.
- Oracle GR-GWR uses the same refit primitive with the true labels.

## Main result

Initial Queen-Ward partition:

- ARI = 0.7671024857
- boundary F1 = 0.2489626556
- estimated boundary edges = 128 vs 113 true

After six accepted refinement iterations:

- ARI = 0.9594008739
- boundary precision = 0.8376068376
- boundary recall = 0.8672566372
- boundary F1 = 0.8521739130
- boundary Jaccard = 0.7424242424
- estimated boundary edges = 117 vs 113 true
- stop reason = no_label_changes

Thus the current refinement substantially improves the pilot segmentation in this controlled single realization.

## Model fit

- OLS: RMSE 0.9534554550, R2 0.0909226953
- GWR: RMSE 0.4551721959, R2 0.7928182720
- MGWR: RMSE 0.4262601675, R2 0.8183022696
- estimated GR-GWR: RMSE 0.4076689414, R2 0.8338060342
- oracle GR-GWR: RMSE 0.4055118406, R2 0.8355601471

Estimated GR-GWR is very close to oracle in in-sample RMSE for this realization, but this is not yet Monte Carlo or out-of-sample evidence.

## Coefficient recovery

Slope RMSE(beta):

| term | GWR | MGWR | estimated GR-GWR | oracle GR-GWR |
|---|---:|---:|---:|---:|
| X1 | 0.228242 | 0.216806 | 0.166985 | 0.037987 |
| X2 | 0.182165 | 0.162769 | 0.137058 | 0.034094 |
| X3 | 0.164201 | 0.145610 | 0.036359 | 0.037482 |

Estimated GR-GWR improves coefficient recovery relative to GWR/MGWR for all three slopes in this realization. X1/X2 still show a meaningful gap from oracle, indicating that remaining boundary-label error matters.

## True-boundary jump recovery

Mean estimated jump / mean true jump across true boundary edges:

| term | GWR | MGWR | estimated GR-GWR | oracle GR-GWR |
|---|---:|---:|---:|---:|
| X1 | 0.144 | 0.165 | 0.760 | 0.953 |
| X2 | 0.178 | 0.236 | 0.852 | 1.060 |
| X3 | 0.173 | 0.255 | 0.971 | 0.984 |

This is the clearest behavioral signal so far: standard GWR/MGWR strongly attenuate the known coefficient discontinuities in this constructed case, whereas estimated GR-GWR recovers most of them and approaches the oracle result.

## Important limitations

- Single realization only.
- K=3 is supplied; unknown-K selection is not tested here.
- The coefficient truth is intentionally favorable to a regime model (piecewise constant within regimes).
- lambda=0.5 is the current working value, not a theoretically selected parameter.
- Reported fit is in-sample.
- No repeated Monte Carlo, no weak-boundary/no-boundary controls, no smooth-within-regime truth, and no irregular topology yet.

## Next simulation tasks

1. Inspect the generated regime and coefficient-recovery figures for pathological behavior.
2. Turn this exact setup into repeated paired Monte Carlo experiments.
3. Add smooth-within / discontinuous-between truth.
4. Add boundary-strength x noise experiments.
5. Add no-boundary negative controls.
6. Add Georgia-topology semi-synthetic experiments.

Result directory:
`results/synthetic/canonical_regime_boundary/`
