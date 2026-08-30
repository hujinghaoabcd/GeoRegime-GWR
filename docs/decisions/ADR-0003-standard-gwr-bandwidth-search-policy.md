# ADR-0003 — Standard GWR bandwidth search policy

Status: **Accepted for the canonical Georgia reproduction baseline**

Date: 2026-08-30

## Context

The research repository needs a self-contained `BasicGWR` that reproduces the canonical Georgia standard-GWR result end to end, including automatic bandwidth selection. The external reference is `mgwr==2.2.1`, but the repository must not depend on `mgwr` for its own fit or automatic search.

While extracting the current PyGWRx adaptive AICc selector, we discovered an important difference in search policy:

- current PyGWRx scans every valid integer adaptive bandwidth and therefore finds the global minimum on the scanned integer AICc curve;
- `mgwr==2.2.1` uses a discrete golden-section search by default.

For the canonical Georgia specification, these policies do not return the same integer:

- current PyGWRx exhaustive search: `k = 116`;
- canonical `mgwr==2.2.1` golden-section search: `k = 117`.

The difference is a search-policy difference, not a failure of the GWR fitting equations. When both implementations are evaluated at `k=117`, local parameters, fitted values, residuals and the hat matrix agree to floating-point machine precision.

## Decision

`BasicGWR(bandwidth="auto")` will use a locally implemented **mgwr-2.2.1-compatible discrete golden-section AICc search** for the canonical standard-GWR benchmark.

The implementation lives in:

- `src/georegime_gwr/bandwidth.py::MGWRCompatibleAICcSelector`

It reproduces the relevant `mgwr==2.2.1` semantics:

- adaptive integer bandwidth;
- initial adaptive search section `40 + 2*p` through `n`, where `p` includes the intercept;
- golden-section constant `0.38197`;
- integer rounding at each evaluation;
- tolerance `1e-6`;
- maximum 200 iterations;
- adaptive compact-kernel boundary distance multiplied by `1.0000001`;
- Gaussian GWR AICc criterion.

The current PyGWRx exhaustive selector is also retained separately as:

- `src/georegime_gwr/bandwidth.py::PyGWRxAdaptiveAICcSelector`

It must not be silently conflated with the benchmark-compatible selector.

## Validation

On the canonical Georgia benchmark:

- external `mgwr.sel_bw.Sel_BW`: bandwidth `117`;
- repository `BasicGWR(bandwidth="auto")`: bandwidth `117`;
- max absolute local-parameter difference: `5.551115123125783e-16`;
- max absolute fitted-value difference: `5.551115123125783e-16`;
- max absolute residual difference: `5.551115123125783e-16`;
- max absolute hat-matrix difference: `5.551115123125783e-17`;
- RSS difference: `-7.105427357601002e-15`;
- AICc difference: `-5.684341886080802e-14`;
- end-to-end validation at tolerance `1e-12`: **PASS**.

Evidence:

- `experiments/validation/basicgwr_vs_mgwr_georgia.py`
- `results/validation/basicgwr_vs_mgwr_georgia/summary.json`
- `.github/workflows/basicgwr-validation.yml`

## Consequences

1. The canonical standard-GWR baseline can be reproduced by repository code without using `mgwr` for its own bandwidth selection or fitting.
2. `mgwr` remains only an external validation oracle in the validation workflow.
3. Search policy is now an explicit methodological choice. A later GR-GWR paper experiment may choose exhaustive AICc search instead, but that change must be documented and must not be presented as identical to the canonical mgwr benchmark.
4. The fact that exhaustive search yields `116` while historical golden-section search yields `117` must be preserved as reproducibility evidence rather than hidden.
