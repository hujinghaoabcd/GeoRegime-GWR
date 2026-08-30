# RELATED_METHODS

> Working bibliography for GR-GWR method positioning. This is a research note, not the final paper literature review.

## Core lineage

### Geographically weighted regression

- Brunsdon, Fotheringham & Charlton (1996), *Geographically Weighted Regression: A Method for Exploring Spatial Nonstationarity*.
- Yu et al. (2020), *On the measurement of bias in geographically weighted regression models*.

Relevance: GWR provides continuous local coefficient surfaces through spatially weighted borrowing. Boundary-crossing borrowing can bias local estimates when neighboring observations belong to different spatial processes.

### Spatial regimes / clustered regression

- Anselin (1990), spatial structural instability / spatial regimes.
- Li & Sang (2019), spatial homogeneity pursuit / spatially clustered coefficients.
- Sugasawa & Murakami (2021), spatially clustered regression.
- Guo, Python & Liu (2023), regionalization algorithms for spatial process heterogeneity.

Relevance: endogenous spatial regimes and spatially contiguous coefficient structures are established ideas; GR-GWR must not claim these as new by themselves.

### Piecewise-smooth / jump-discontinuous coefficient fields

- Zhu, Fan & Kong (2014), spatially varying coefficient models with jump discontinuities.
- Lin et al. (2022), *Spatially Clustered Varying Coefficient Model*.

Relevance: “smooth within regions, discontinuous across boundaries” is already established statistically and cannot be claimed as GR-GWR's original theoretical idea.

### GWR coefficient clustering

- Lee, Gangnon & Zhu (2017), cluster detection of spatial regression coefficients.
- Nicholson et al. (2019), GWR followed by spatial contiguity-constrained clustering.
- Areed et al. (2024), Bayesian clustered coefficients / Bayesian cluster GWR work.

Relevance: `GWR -> local coefficients -> spatial clustering` already exists. GR-GWR's potential contribution must lie beyond post-hoc clustering.

### Boundary-aware GWR

- Sheng et al. (2026), multiscale GWR considering boundary effects.

Relevance: boundary-aware spatial borrowing is already an active GWR extension. A key distinction to investigate is exogenous/known boundaries versus endogenous regime discovery.

## Current novelty hypothesis

The current working novelty hypothesis is not any single component above. It is the feedback coupling:

`local GWR relationships -> endogenous contiguous regimes -> regime-gated local GWR -> predictive regime refinement`

This remains a hypothesis until systematic literature checking and experiments are complete.

## Literature questions still open

1. Is there an earlier method with the same full feedback loop?
2. Are there GWR methods that directly accept polygon contiguity / spatial weights `W` for regime discovery?
3. Which existing methods are realistically reproducible competitors for the paper experiments?
4. How should statistical inference account for data-driven regime selection?
