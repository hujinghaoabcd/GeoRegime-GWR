"""Plot the canonical Georgia county polygons and the Queen-contiguity graph.

This script is intentionally an experiment/visualization utility, not a package API.
It consumes the same shapefile and the persisted Queen matrix used by the GR-GWR
spatial-topology experiment so the plotted graph is exactly the W matrix that will
be used in later algorithm steps.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
W_CSV = ROOT / "results" / "spatial" / "georgia_queen_w" / "queen_adjacency_matrix.csv"
OUT = ROOT / "results" / "spatial" / "georgia_queen_w"


def _find_column(columns, target: str):
    target = target.lower()
    for col in columns:
        if str(col).lower() == target:
            return col
    return None


def _load_inputs():
    gdf = gpd.read_file(SHP)
    if len(gdf) != 159:
        raise RuntimeError(f"Expected 159 Georgia counties, found {len(gdf)}")

    area_col = _find_column(gdf.columns, "AreaKey")
    if area_col is None:
        raise RuntimeError("AreaKey column not found in Georgia shapefile")

    area_keys = gdf[area_col].astype(str).tolist()
    W_df = pd.read_csv(W_CSV, index_col=0)
    W_df.index = W_df.index.astype(str)
    W_df.columns = W_df.columns.astype(str)

    if area_keys != W_df.index.tolist() or area_keys != W_df.columns.tolist():
        raise RuntimeError("Shapefile row order and Queen W AreaKey order do not match")

    W = W_df.to_numpy(dtype=int)
    if W.shape != (159, 159):
        raise RuntimeError(f"Unexpected W shape: {W.shape}")
    if not np.array_equal(W, W.T):
        raise RuntimeError("Queen W must be symmetric")

    return gdf, W


def _draw_counties(ax, gdf):
    gdf.plot(ax=ax, facecolor="white", edgecolor="black", linewidth=0.55)
    ax.set_aspect("equal")
    ax.set_axis_off()


def _draw_graph(ax, gdf, W):
    _draw_counties(ax, gdf)
    centers = gdf.geometry.centroid

    for i in range(len(gdf)):
        xi, yi = centers.iloc[i].x, centers.iloc[i].y
        for j in range(i + 1, len(gdf)):
            if W[i, j] == 1:
                xj, yj = centers.iloc[j].x, centers.iloc[j].y
                ax.plot([xi, xj], [yi, yj], linewidth=0.42, alpha=0.42, zorder=2)

    ax.scatter(centers.x, centers.y, s=7, zorder=3)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gdf, W = _load_inputs()
    edge_count = int(W.sum() // 2)

    # 1. Pure county map.
    fig, ax = plt.subplots(figsize=(8.5, 9.0))
    _draw_counties(ax, gdf)
    ax.set_title("Georgia counties", fontsize=15)
    fig.tight_layout()
    fig.savefig(OUT / "georgia_counties_map.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "georgia_counties_map.svg", bbox_inches="tight")
    plt.close(fig)

    # 2. Queen-contiguity graph over the polygons.
    fig, ax = plt.subplots(figsize=(8.5, 9.0))
    _draw_graph(ax, gdf, W)
    ax.set_title(f"Georgia Queen-contiguity graph (159 counties, {edge_count} edges)", fontsize=14)
    fig.tight_layout()
    fig.savefig(OUT / "georgia_queen_graph.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "georgia_queen_graph.svg", bbox_inches="tight")
    plt.close(fig)

    # 3. Side-by-side figure for direct visual inspection.
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 8.2))
    _draw_counties(axes[0], gdf)
    axes[0].set_title("(a) County polygons", fontsize=14)

    _draw_graph(axes[1], gdf, W)
    axes[1].set_title("(b) Queen adjacency W", fontsize=14)

    fig.suptitle("Georgia spatial structure used by GR-GWR", fontsize=16)
    fig.tight_layout()
    fig.savefig(OUT / "georgia_spatial_structure.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "georgia_spatial_structure.svg", bbox_inches="tight")
    plt.close(fig)

    print(f"Generated Georgia spatial maps in: {OUT.relative_to(ROOT)}")
    print(f"Counties: {len(gdf)}")
    print(f"Queen undirected edges: {edge_count}")


if __name__ == "__main__":
    main()
