"""Build a binary symmetric Queen-contiguity matrix for the Georgia counties.

This is a reproducible research-data construction step for GR-GWR. It does
not define package APIs. The output matrix is ordered by the shapefile row
order and labelled with AreaKey values so it can be joined back to the
canonical Georgia regression table.
"""

from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from libpysal.weights import Queen
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[2]
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
OUT = ROOT / "results" / "spatial" / "georgia_queen_w"


def _find_column(columns, target: str):
    target = target.lower()
    for col in columns:
        if str(col).lower() == target:
            return col
    return None


def main() -> None:
    gdf = gpd.read_file(SHP)
    n = len(gdf)
    if n != 159:
        raise RuntimeError(f"Expected 159 Georgia counties, found {n}")

    area_col = _find_column(gdf.columns, "AreaKey")
    if area_col is None:
        raise RuntimeError(f"AreaKey column not found; columns={list(gdf.columns)}")
    name_col = _find_column(gdf.columns, "Name")

    area_keys = gdf[area_col].astype(str).tolist()
    if len(set(area_keys)) != n:
        raise RuntimeError("AreaKey must be unique")

    queen = Queen.from_dataframe(gdf, use_index=False)

    W = np.zeros((n, n), dtype=np.uint8)
    for i, neigh in queen.neighbors.items():
        for j in neigh:
            W[int(i), int(j)] = 1

    W = np.maximum(W, W.T)
    np.fill_diagonal(W, 0)

    if not np.array_equal(W, W.T):
        raise RuntimeError("Queen W is not symmetric")
    if np.any(np.diag(W) != 0):
        raise RuntimeError("Queen W has self-neighbours")
    if not np.all(np.isin(W, [0, 1])):
        raise RuntimeError("Queen W is not binary")

    degrees = W.sum(axis=1).astype(int)
    edges = int(W.sum() // 2)
    n_components, labels = connected_components(
        csr_matrix(W), directed=False, return_labels=True
    )
    component_sizes = np.bincount(labels).astype(int).tolist()
    isolates = np.flatnonzero(degrees == 0).tolist()

    OUT.mkdir(parents=True, exist_ok=True)

    matrix = pd.DataFrame(W, index=area_keys, columns=area_keys)
    matrix.index.name = "AreaKey"
    matrix.to_csv(OUT / "queen_adjacency_matrix.csv")

    degree_df = pd.DataFrame({
        "row_index": np.arange(n, dtype=int),
        "AreaKey": area_keys,
        "degree": degrees,
        "component": labels.astype(int),
    })
    if name_col is not None:
        degree_df.insert(2, "Name", gdf[name_col].astype(str).to_numpy())
    degree_df.to_csv(OUT / "degree_summary.csv", index=False)

    edge_rows = []
    for i in range(n):
        for j in range(i + 1, n):
            if W[i, j]:
                row = {
                    "i": i,
                    "j": j,
                    "AreaKey_i": area_keys[i],
                    "AreaKey_j": area_keys[j],
                }
                if name_col is not None:
                    row["Name_i"] = str(gdf.iloc[i][name_col])
                    row["Name_j"] = str(gdf.iloc[j][name_col])
                edge_rows.append(row)
    pd.DataFrame(edge_rows).to_csv(OUT / "edge_list.csv", index=False)

    preview_n = 20
    matrix.iloc[:preview_n, :preview_n].to_csv(OUT / "matrix_preview_20.csv")

    summary = {
        "source_shapefile": str(SHP.relative_to(ROOT)),
        "construction": "Queen contiguity",
        "n": n,
        "matrix_shape": [n, n],
        "binary": True,
        "symmetric": True,
        "zero_diagonal": True,
        "undirected_edge_count": edges,
        "possible_undirected_edges": int(n * (n - 1) // 2),
        "density": float(edges / (n * (n - 1) / 2)),
        "degree_min": int(degrees.min()),
        "degree_max": int(degrees.max()),
        "degree_mean": float(degrees.mean()),
        "degree_median": float(np.median(degrees)),
        "isolated_count": len(isolates),
        "isolated_AreaKeys": [area_keys[i] for i in isolates],
        "connected_components": int(n_components),
        "component_sizes": component_sizes,
        "AreaKey_order": area_keys,
    }
    (OUT / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("\nFirst 20x20 Queen adjacency block:\n")
    print(matrix.iloc[:preview_n, :preview_n].to_string())


if __name__ == "__main__":
    main()
