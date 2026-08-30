"""Fit the repository's own BasicGWR on Georgia and map local coefficients.

This script deliberately uses ``georegime_gwr.gwr.BasicGWR`` only for model
fitting.  It does not import or run mgwr.GWR.

The regression is the same standardized Georgia baseline used by the research
validation path:

    PctBach ~ PctFB + PctBlack + PctRural

The research-default adaptive exhaustive AICc selector is expected to choose
k=116.  Outputs include one 2x2 figure (PNG + SVG) containing the intercept and
three local slope surfaces, plus county-level and summary CSV tables.
"""

from __future__ import annotations

from pathlib import Path
import sys

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr.gwr import BasicGWR  # noqa: E402

CSV = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
OUT = ROOT / "results" / "real_data" / "georgia_gwr_baseline"

TERMS = ["Intercept", "PctFB", "PctBlack", "PctRural"]
MAP_COLUMNS = {
    "Intercept": "coef_Intercept",
    "PctFB": "coef_PctFB",
    "PctBlack": "coef_PctBlack",
    "PctRural": "coef_PctRural",
}
TITLES = {
    "Intercept": "(a) Intercept",
    "PctFB": "(b) PctFB coefficient",
    "PctBlack": "(c) PctBlack coefficient",
    "PctRural": "(d) PctRural coefficient",
}


def _find_column(columns, target: str):
    target = target.lower()
    for col in columns:
        if str(col).lower() == target:
            return col
    return None


def _fit_basicgwr(df: pd.DataFrame) -> BasicGWR:
    y = df["PctBach"].to_numpy(dtype=float).reshape(-1, 1)
    X = df[["PctFB", "PctBlack", "PctRural"]].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)

    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean(axis=0)) / y.std(axis=0, ddof=0)

    model = BasicGWR(
        bandwidth="auto",
        kernel="bisquare",
        fit_intercept=True,
    ).fit(Xz, yz.reshape(-1), coords)

    if int(model.bandwidth_) != 116:
        raise RuntimeError(
            f"Expected research-default Georgia bandwidth k=116; got {model.bandwidth_}"
        )
    return model


def _attach_coefficients(
    gdf: gpd.GeoDataFrame,
    df: pd.DataFrame,
    params: np.ndarray,
) -> gpd.GeoDataFrame:
    shp_key = _find_column(gdf.columns, "AreaKey")
    csv_key = _find_column(df.columns, "AreaKey")
    if shp_key is None or csv_key is None:
        raise RuntimeError("AreaKey column not found in Georgia shapefile or CSV")

    coef = pd.DataFrame({"AreaKey_join": df[csv_key].astype(str)})
    for j, term in enumerate(TERMS):
        coef[MAP_COLUMNS[term]] = params[:, j]

    out = gdf.copy()
    out["AreaKey_join"] = out[shp_key].astype(str)
    out = out.merge(coef, on="AreaKey_join", how="left", validate="1:1", sort=False)

    missing = list(MAP_COLUMNS.values())
    if out[missing].isna().any().any():
        raise RuntimeError("Missing local coefficients after joining them to Georgia polygons")
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(CSV)
    gdf = gpd.read_file(SHP)
    if len(df) != 159 or len(gdf) != 159:
        raise RuntimeError(
            f"Expected 159 Georgia counties; csv={len(df)}, shapefile={len(gdf)}"
        )

    model = _fit_basicgwr(df)
    params = np.asarray(model.parameters_, dtype=float)
    if params.shape != (159, 4):
        raise RuntimeError(f"Expected coefficient matrix shape (159, 4); got {params.shape}")

    mapped = _attach_coefficients(gdf, df, params)

    # Persist county-level coefficients in regression-data order.
    coef_table = pd.DataFrame({"AreaKey": df["AreaKey"]})
    for j, term in enumerate(TERMS):
        coef_table[term] = params[:, j]
    coef_table.to_csv(OUT / "basicgwr_local_coefficients.csv", index=False)

    summary = pd.DataFrame(
        {
            "term": TERMS,
            "min": [float(np.min(params[:, j])) for j in range(4)],
            "mean": [float(np.mean(params[:, j])) for j in range(4)],
            "std": [float(np.std(params[:, j], ddof=0)) for j in range(4)],
            "max": [float(np.max(params[:, j])) for j in range(4)],
        }
    )
    summary.to_csv(OUT / "basicgwr_local_coefficients_summary.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 13.2))

    for ax, term in zip(axes.flat, TERMS):
        col = MAP_COLUMNS[term]
        values = mapped[col].to_numpy(dtype=float)
        limit = float(np.max(np.abs(values)))
        if limit <= 0.0:
            limit = 1.0
        norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)

        mapped.plot(
            column=col,
            ax=ax,
            cmap="coolwarm",
            norm=norm,
            legend=True,
            edgecolor="black",
            linewidth=0.28,
            legend_kwds={"shrink": 0.72},
        )
        ax.set_title(TITLES[term], fontsize=14)
        ax.set_aspect("equal")
        ax.set_axis_off()

    fig.suptitle(
        "Georgia BasicGWR local coefficients (research default, adaptive k=116)",
        fontsize=17,
        y=0.995,
    )
    fig.tight_layout()

    png = OUT / "georgia_basicgwr_local_coefficients.png"
    svg = OUT / "georgia_basicgwr_local_coefficients.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print(f"BasicGWR bandwidth: {int(model.bandwidth_)}")
    print(summary.to_string(index=False))
    print(f"Generated: {png.relative_to(ROOT)}")
    print(f"Generated: {svg.relative_to(ROOT)}")
    print(f"Generated: {(OUT / 'basicgwr_local_coefficients.csv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
