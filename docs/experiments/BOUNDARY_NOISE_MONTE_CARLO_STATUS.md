# Paired Boundary-Noise Monte Carlo Status

Last updated: 2026-09-01

## Design

Publication-oriented paired Monte Carlo follow-up to Simulation 3.

- lattice: 25 x 25 (n=625)
- Queen topology
- current estimated chain: pilot GWR -> Queen-constrained Ward -> lambda=0.5 boundary refinement -> unified RegimeAwareGWR
- K=3 fixed in these experiments to isolate boundary detectability from K selection
- four models use exactly the same generated data at each seed: GWR, MGWR, estimated GR-GWR, Oracle GR-GWR
- 5 representative scenarios x 100 independent seeds = 500 simulated data sets
- all 500 runs completed successfully; MGWR also succeeded in all 500 runs
- primary evidence: coefficient RMSE against known truth, paired 95% confidence intervals, recovery rates; in-sample fit is secondary only

Scenarios:

1. `null_forced_k3`: delta=0, noise=0.75
2. `weak`: delta=0.5, noise=0.75
3. `moderate`: delta=1.0, noise=0.75
4. `strong`: delta=2.0, noise=0.35
5. `high_noise_failure`: delta=1.5, noise=2.0 (name inherited from the single-seed prototype; Monte Carlo shows it is not a stable failure case)

## Main coefficient-recovery results

### Strong boundary

Mean slope coefficient RMSE:

- GWR: 0.17095
- MGWR: 0.15283
- GR-GWR: 0.03786
- Oracle: 0.01229

GR-GWR reduction relative to GWR: 77.85%.
GR-GWR reduction relative to MGWR: 75.22%.

Paired GR-GWR minus GWR coefficient-RMSE difference: -0.13309, 95% CI [-0.13905, -0.12713].
Paired GR-GWR minus MGWR: -0.11496, 95% CI [-0.12085, -0.10908].
GR-GWR coefficient-RMSE win rate: 100% vs GWR and 100% vs MGWR.

Boundary recovery:

- refined ARI mean: 0.99212
- refined boundary F1 mean: 0.95863
- GR-GWR true-jump recovery mean: 0.96519
- Oracle jump recovery mean: 1.00185

Interpretation: when discontinuities are strong and identifiable, estimated GR-GWR nearly recovers the true regime structure and coefficient jumps and is close to the oracle boundary-aware fit.

### Moderate boundary

Mean slope coefficient RMSE:

- GWR: 0.16428
- MGWR: 0.16051
- GR-GWR: 0.10936
- Oracle: 0.04091

GR-GWR reduction: 33.43% vs GWR and 31.87% vs MGWR.
Paired differences are strictly negative at 95% confidence:

- vs GWR: -0.05493, 95% CI [-0.06095, -0.04890]
- vs MGWR: -0.05115, 95% CI [-0.05719, -0.04511]

Coefficient-RMSE win rates: 96% vs GWR, 94% vs MGWR.
Boundary recovery: refined ARI 0.93028; boundary F1 0.70727; jump recovery 0.74301.

### Weak boundary

Mean slope coefficient RMSE:

- GWR: 0.13104
- MGWR: 0.13168
- GR-GWR: 0.12293
- Oracle: 0.05367

GR-GWR reduction: 6.19% vs GWR and 6.64% vs MGWR.
Although the gain is much smaller, paired confidence intervals still exclude zero:

- vs GWR: -0.00811, 95% CI [-0.01197, -0.00425]
- vs MGWR: -0.00875, 95% CI [-0.01259, -0.00490]

Coefficient-RMSE win rates: 68% vs GWR, 69% vs MGWR.
Boundary recovery is substantially weaker: refined ARI 0.80628; boundary F1 0.40135; jump recovery 0.45586.

Interpretation: weak discontinuities remain detectable on average but performance becomes seed-sensitive and the gain over smooth GWR/MGWR is modest.

### High-noise scenario

The single-seed phase diagram had suggested this could be a failure case. The 100-seed paired Monte Carlo does **not** support a stable failure conclusion.

Mean slope coefficient RMSE:

- GWR: 0.16109
- MGWR: 0.16153
- GR-GWR: 0.14847
- Oracle: 0.06063

GR-GWR reduction: 7.83% vs GWR and 8.09% vs MGWR.
Paired differences:

- vs GWR: -0.01262, 95% CI [-0.01788, -0.00736]
- vs MGWR: -0.01306, 95% CI [-0.01829, -0.00784]

Coefficient-RMSE win rates: 68% vs GWR, 69% vs MGWR.
Boundary recovery remains mediocre: refined ARI 0.83613; boundary F1 0.43405; jump recovery 0.49434.

Interpretation: high noise weakens boundary identification and leaves a large gap to the oracle, but under this specific delta=1.5 setting estimated GR-GWR still improves coefficient recovery on average. The previous single-seed failure was therefore not robust.

## Critical null result: no discontinuity but K=3 forced

This is the strongest negative-control result.

When delta=0 there is no coefficient discontinuity, but the prototype is forced to fit K=3.

Mean slope coefficient RMSE:

- GWR: 0.03529
- MGWR: 0.04577
- forced-K GR-GWR: 0.07802
- Oracle K=3 boundary-aware fit: 0.06344

Thus forced GR-GWR coefficient RMSE is 121.1% higher than GWR and 70.5% higher than MGWR.
Paired differences are decisively positive:

- GR-GWR minus GWR: +0.04274, 95% CI [+0.04006, +0.04542]
- GR-GWR minus MGWR: +0.03225, 95% CI [+0.02902, +0.03549]

Coefficient-RMSE win rate: 0% vs GWR and only 3% vs MGWR.

At the same time, forced GR-GWR improves in-sample fitted RMSE in 100% of runs versus GWR and 93% versus MGWR.

This is direct Monte Carlo evidence that training fit can reward false segmentation. A publication-ready GR-GWR must therefore allow K=1 / no-boundary selection and must not choose segmentation using in-sample fit alone.

## In-sample fit results (secondary diagnostic)

GR-GWR mean fitted RMSE is lower than GWR and MGWR in all five scenarios, including the null forced-K case. This confirms that fitted RMSE alone is not a valid segmentation criterion. Known-truth coefficient recovery and ultimately spatial out-of-sample validation are required.

## Runtime diagnostic (n=625)

Mean runtime per replication:

- GWR: approximately 1.39-1.54 s across scenarios
- MGWR: approximately 14.32-31.53 s
- full GR-GWR chain including pilot GWR: approximately 1.86-2.19 s
- GR refinement/refit after the pilot: approximately 0.45-0.65 s

This implementation examines current boundary candidates rather than all spatial units at every refinement step. The result is encouraging but is not a formal scaling proof; computational scaling still requires a dedicated n-sensitivity experiment.

## Current publication-level conclusions

1. The strong-boundary advantage is highly reproducible across random seeds.
2. The moderate-boundary advantage is also strong and statistically stable.
3. Weak boundaries yield smaller and more seed-sensitive gains, but the mean paired coefficient-recovery advantage remains significant in this design.
4. High noise reduces boundary accuracy and the gap to Oracle remains substantial; nevertheless the earlier single-seed failure does not replicate as a systematic failure.
5. The real failure mode demonstrated robustly is false segmentation when no true discontinuity exists.
6. Therefore K=1 / no-boundary selection is not optional; it is a core requirement of the final method.
7. In-sample fit cannot select K or justify segmentation.
8. Boundary identification quality is strongly linked to the magnitude of the GR-GWR advantage.

## Next methodological task

Do not add more forced-K Monte Carlo grids yet. The next highest-priority experiment is the formal no-boundary / K-selection control in which candidate K includes K=1. The method must demonstrate that stationary or purely smooth truth can revert to ordinary GWR instead of inventing regimes. After that, test K recovery when K_true > 1 and then proceed to out-of-sample spatial validation.