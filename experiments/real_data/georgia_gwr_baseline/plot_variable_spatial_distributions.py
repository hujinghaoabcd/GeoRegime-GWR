"""Plot Georgia baseline GWR variables as one multi-panel county map.

The figure uses the canonical Georgia polygon layer and the exact regression
variables used by the standard-GWR baseline experiment:

- PctBach  : dependent variable
- PctFB    : percentage foreign born
- PctBlack : percentage Black population
- PctRural : percentage rural population

This is an experiment/visualization utility, not a package API.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
CSV = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
OUT = ROOT / "results" / "real_data" / "georgia_gwr_baseline"

VARIABLES = ["PctBach", "PctFB", "PctBlack", "PctRural"]
TITLES = {
    "PctBach": "(a) PctBach",
    "PctFB": "(b) PctFB",
    "PctBlack": "(c) PctBlack",
    "PctRural": "(d) PctRural",
}


def _find_column(columns, target: str):
    target = target.lower()
    for col in columns:
        if str(col).lower() == target:
            return col
    return None


def _load_data() -> gpd.GeoDataFrame:
    gdf = gpd.read_file(SHP)
    df = pd.read_csv(CSV)

    if len(gdf) != 159 or len(df) != 159:
        raise RuntimeError(
            f"Expected 159 Georgia counties; shapefile={len(gdf)}, csv={len(df)}"
        )

    shp_key = _find_column(gdf.columns, "AreaKey")
    csv_key = _find_column(df.columns, "AreaKey")
    if shp_key is None or csv_key is None:
        raise RuntimeError("AreaKey column not found in Georgia shapefile or CSV")

    gdf = gdf.copy()
    df = df.copy()
    gdf[shp_key] = gdf[shp_key].astype(str)
    df[csv_key] = df[csv_key].astype(str)

    # The canonical Georgia shapefile already contains the baseline variables.
    # If a variable is ever absent there, fill only that missing variable from
    # GData_utm.csv by AreaKey. This avoids pandas _x/_y merge suffixes.
    csv_by_key = df.set_index(csv_key)
    for variable in VARIABLES:
        if variable not in gdf.columns:
            if variable not in csv_by_key.columns:
                raise RuntimeError(f"Missing baseline variable: {variable}")
            gdf[variable] = gdf[shp_key].map(csv_by_key[variable])

        gdf[variable] = pd.to_numeric(gdf[variable], errors="coerce")

    if gdf[VARIABLES].isna().any().any():
        raise RuntimeError("Missing or non-numeric baseline values in Georgia polygons")

    return gdf


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gdf = _load_data()

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 14.0))

    for ax, variable in zip(axes.flat, VARIABLES):
        gdf.plot(
            column=variable,
            ax=ax,
            legend=True,
            edgecolor="black",
            linewidth=0.28,
            legend_kwds={"shrink": 0.72},
        )
        ax.set_title(TITLES[variable], fontsize=14)
        ax.set_aspect("equal")
        ax.set_axis_off()

    fig.suptitle(
        "Georgia county-level spatial distributions of baseline GWR variables",
        fontsize=17,
        y=0.995,
    )
    fig.tight_layout()

    png = OUT / "georgia_variable_spatial_distributions.png"
    svg = OUT / "georgia_variable_spatial_distributions.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    print(f"Generated: {png.relative_to(ROOT)}")
    print(f"Generated: {svg.relative_to(ROOT)}")
    print(f"Variables: {', '.join(VARIABLES)}")


if __name__ == "__main__":
    main()
