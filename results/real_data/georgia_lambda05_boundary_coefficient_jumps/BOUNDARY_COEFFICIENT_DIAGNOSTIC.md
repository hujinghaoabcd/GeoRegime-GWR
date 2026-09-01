# Lambda=0.5 boundary coefficient diagnostic

Status: exploratory/descriptive only. The final regimes are data-adaptive and learned from the same response data, so these diagnostics are not independent validation of true boundaries.

## Setup

- Georgia n = 159 counties.
- Final working GR-GWR partition at lambda=0.5: sizes 62, 14, 12, 38, 16, 17.
- Queen graph: 431 undirected edges.
- Final cross-regime boundary edges: 71.
- Within-regime Queen edges: 360.
- Models compared: ordinary GWR, external MGWR, working GR-GWR.
- Slopes: PctFB, PctBlack, PctRural.

## Main result: cross-boundary coefficient jumps

Mean absolute jump across the 71 final boundary edges:

| Term | GWR | MGWR | GR-GWR | GR/GWR |
|---|---:|---:|---:|---:|
| PctFB | 0.0616 | 0.0517 | 0.4271 | 6.93x |
| PctBlack | 0.0302 | 0.0134 | 0.1300 | 4.31x |
| PctRural | 0.0219 | 0.0054 | 0.3328 | 15.17x |

Fraction of final boundary edges where GR-GWR has a larger coefficient jump than ordinary GWR:

- PctFB: 97.18%.
- PctBlack: 85.92%.
- PctRural: 97.18%.

## Boundary-vs-interior contrast

Mean absolute jump ratio, boundary edges divided by within-regime edges:

| Term | GWR | MGWR | GR-GWR |
|---|---:|---:|---:|
| PctFB | 2.07 | 1.82 | 8.26 |
| PctBlack | 1.61 | 1.44 | 5.98 |
| PctRural | 1.82 | 1.39 | 5.95 |

Thus the increase in local coefficient variation under GR-GWR is strongly concentrated at inferred regime boundaries rather than being a uniform amplification everywhere. However, within-regime jumps also increase for some terms, especially PctRural, so the result must not be described as pure boundary recovery.

## Interpretation

The diagnostic is consistent with the intended GR-GWR mechanism: ordinary GWR/MGWR retain smooth distance-based coefficient surfaces across the inferred boundaries, whereas the regime-aware fit permits stronger local coefficient discontinuities where cross-regime borrowing is prohibited.

This is not yet proof that the inferred boundaries are true. The partition itself was learned adaptively from the same response data, and the hard regime mask mechanically encourages discontinuity across the selected boundary. Independent synthetic experiments with known true boundaries and spatial/out-of-sample validation are still required.

Generated figures:

- `PctFB_gwr_grgwr_delta_map.png/.svg`
- `PctBlack_gwr_grgwr_delta_map.png/.svg`
- `PctRural_gwr_grgwr_delta_map.png/.svg`
- `cross_boundary_jump_ratio.png/.svg`
