# ADR-0003 — Standard GWR bandwidth search policy

Status: **Accepted**

Date: 2026-08-30

## Context

The research repository needs a self-contained `BasicGWR` for two distinct purposes:

1. provide a rigorous standard-GWR baseline for later GR-GWR paper experiments;
2. reproduce the canonical `mgwr==2.2.1` Georgia standard-GWR example when historical compatibility is explicitly requested.

These purposes must not be conflated.

Current PyGWRx and historical `mgwr==2.2.1` use different adaptive-bandwidth search policies:

- current PyGWRx treats adaptive bandwidth as a discrete integer neighbour order and exhaustively evaluates every valid `k`;
- `mgwr==2.2.1` uses a rounded discrete golden-section search by default.

On canonical Georgia these policies differ:

- exhaustive integer AICc minimum: `k = 116`, AICc about `298.9856`;
- historical mgwr discrete golden-section result: `k = 117`, AICc about `299.0508`.

When both fit GWR at the same `k=117`, the local coefficients, fitted values, residuals and hat matrix agree to floating-point machine precision. The 116/117 difference is therefore a search-policy difference, not a GWR-equation discrepancy.

## Decision

The repository exposes **three explicit search modes**.

### 1. Adaptive research default — exhaustive integer AICc

`BasicGWR(bandwidth="auto")`

with adaptive bandwidth (the default for `"auto"`) uses:

- `PyGWRxAdaptiveAICcSelector`;
- every valid integer `k` in the search domain is evaluated;
- the minimum finite AICc is returned;
- this is the **default research policy** because it gives the global minimum on the specified discrete candidate domain.

For Georgia this returns `k=116`.

### 2. Fixed-distance research default — golden-section AICc

`BasicGWR(bandwidth="auto", adaptive=False)`

uses:

- `FixedGoldenAICcSelector`;
- continuous golden-section minimization over a positive distance interval;
- this reflects the fact that a fixed bandwidth is continuous rather than an integer neighbour order.

The fixed-distance selector is implemented, but it has not yet been externally benchmarked against a dedicated fixed-GWR reference fixture. It must not be described as externally validated until such a test is added.

### 3. Historical mgwr reproduction — explicit compatibility mode

`BasicGWR(
    bandwidth="auto",
    adaptive=True,
    search_strategy="mgwr_golden",
)`

uses:

- `MGWRCompatibleAICcSelector`;
- `mgwr==2.2.1` adaptive initial section semantics;
- golden-section constant `0.38197`;
- integer rounding at candidate evaluations;
- tolerance `1e-6` and maximum 200 iterations;
- historical compact-kernel boundary inflation `* 1.0000001`;
- Gaussian GWR AICc.

For canonical Georgia this mode is expected to return `k=117` and is used only when exact historical/canonical reproduction is the goal.

## Rationale

Adaptive bandwidth is an integer decision variable:

`k ∈ {k_min, ..., n}`.

For paper experiments, exhaustive evaluation removes optimizer uncertainty from the standard-GWR baseline. A reviewer can therefore not attribute a GR-GWR improvement to a standard GWR baseline whose discrete bandwidth optimizer stopped at a non-global candidate.

Fixed-distance bandwidth is continuous, so exhaustive enumeration is not meaningful; golden-section or another continuous one-dimensional optimizer is appropriate.

The historical mgwr search remains valuable for reproducibility, but reproducibility mode must not silently redefine the research default.

## Implementation

- `src/georegime_gwr/bandwidth.py::PyGWRxAdaptiveAICcSelector`
- `src/georegime_gwr/bandwidth.py::FixedGoldenAICcSelector`
- `src/georegime_gwr/bandwidth.py::MGWRCompatibleAICcSelector`
- `src/georegime_gwr/gwr.py::BasicGWR`

The selected policy and search result are exposed through fitted attributes such as `search_strategy_`, `bandwidth_`, and `bandwidth_search_`.

## Validation contract

The canonical Georgia validation must simultaneously check that:

1. external `mgwr.sel_bw.Sel_BW` returns `117`;
2. repository adaptive research default returns `116`;
3. explicit `mgwr_golden` compatibility mode returns `117`;
4. compatibility-mode GWR at 117 matches external `mgwr.gwr.GWR` numerically to machine precision;
5. the exhaustive 116 solution has lower AICc than the compatibility 117 solution.

Evidence is maintained in:

- `experiments/validation/basicgwr_vs_mgwr_georgia.py`
- `results/validation/basicgwr_vs_mgwr_georgia/`
- `.github/workflows/basicgwr-validation.yml`

## Consequences

1. **Do not say that `BasicGWR(bandwidth="auto")` reproduces mgwr 117 by default.** The adaptive research default is exhaustive and returns 116 on Georgia.
2. Use `search_strategy="mgwr_golden"` only for historical/canonical mgwr reproduction.
3. Standard-GWR paper comparisons should normally use the exhaustive adaptive policy unless a later ADR explicitly changes the experimental protocol.
4. The 116/117 distinction is reproducibility evidence and must remain documented.
5. GR-GWR's final bandwidth policy remains a separate methodological decision and is not frozen by this ADR.
