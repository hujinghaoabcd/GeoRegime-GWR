# Georgia GWR/MGWR reproduction

This experiment reproduces the canonical Georgia example used by **Oshan et al. (2019), _mgwr: A Python Implementation of Multiscale Geographically Weighted Regression for Investigating Process Spatial Heterogeneity and Scale_**.

## Data

Canonical PySAL/libpysal Georgia data:

- 159 Georgia counties;
- 1990 U.S. Census socio-demographic variables;
- projected county coordinates `X`, `Y`;
- county polygons `G_utm.*` retained for later spatial-weights and mapping work.

The source copy is vendored into `data/raw/georgia/` by the reproduction workflow.

## Model specification

Dependent variable:

- `PctBach`: percentage of residents with a bachelor's degree or higher.

Explanatory variables, in this exact order:

1. `PctFB`: percentage foreign-born;
2. `PctBlack`: percentage identifying as African American;
3. `PctRural`: percentage classified as rural.

Both `y` and the three X columns are standardized with population standard deviation (`ddof=0`), matching the published PySAL notebook.

### GWR

- Gaussian model;
- adaptive bisquare kernel;
- AICc bandwidth search;
- historical notebook bandwidth: **117**.

### MGWR

- Gaussian model;
- adaptive bisquare kernels;
- covariate-specific AICc bandwidth search;
- historical notebook bandwidths for `[Intercept, PctFB, PctBlack, PctRural]`: **[92, 101, 136, 158]**.

## Why this experiment exists

This is a **trusted external baseline**, not a GR-GWR experiment yet. It serves three purposes:

1. confirm that the canonical Georgia data are understood correctly;
2. establish reference GWR/MGWR outputs before modifying GR-GWR;
3. provide a real-data benchmark that later GR-GWR variants can use without changing the data/model specification silently.

## Run

```bash
python -m pip install "mgwr==2.2.1" "libpysal>=4.10,<5" pandas
python experiments/real_data/georgia_mgwr_2019/reproduce.py --vendor-data
```

Outputs are written to `results/reproduction/georgia_mgwr_2019/`.
