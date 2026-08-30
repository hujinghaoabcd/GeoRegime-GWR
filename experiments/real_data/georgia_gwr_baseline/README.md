# Georgia standard-GWR baseline

当前阶段只研究 **标准 GWR**。本实验虽然使用 `mgwr` Python 包，但只调用其中的 `mgwr.gwr.GWR`，**不运行 MGWR**。

## Data

Canonical PySAL/libpysal Georgia data:

- 159 Georgia counties;
- 1990 U.S. Census socio-demographic variables;
- projected county coordinates `X`, `Y`;
- county polygons `G_utm.*` retained for later Queen/Rook adjacency experiments.

The canonical copy is stored under `data/raw/georgia/`.

## Model specification

Dependent variable:

- `PctBach`: percentage with a bachelor's degree or higher.

Explanatory variables, in this exact order:

1. `PctFB`: percentage foreign-born;
2. `PctBlack`: percentage Black;
3. `PctRural`: percentage rural.

Pre-processing:

- z-score `PctBach`;
- z-score each explanatory variable;
- population standard deviation (`ddof=0`).

GWR settings:

- Gaussian response;
- adaptive bisquare kernel;
- intercept included;
- AICc bandwidth selection;
- canonical reference bandwidth: **117 neighbours**.

## Purpose

This experiment is the trusted **external standard-GWR baseline** for the GR-GWR paper repository.

Research order:

1. reproduce canonical GWR with `mgwr.gwr.GWR`;
2. compare repository `BasicGWR` against this result point-by-point;
3. only after the standard GWR baseline is validated, modify or test GR-GWR.

MGWR and multiscale bandwidths are deliberately out of scope for this phase.

## Run

```bash
python -m pip install "mgwr==2.2.1" "libpysal==4.10.0" "numpy==1.26.4" "scipy==1.11.4" "pandas==2.2.3"
python experiments/real_data/georgia_gwr_baseline/reproduce.py --vendor-data
```

Outputs are written to:

`results/reproduction/georgia_gwr_baseline/`
