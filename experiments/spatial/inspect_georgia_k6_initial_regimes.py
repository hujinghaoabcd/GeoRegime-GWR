"""Focused inspection of the exploratory Georgia K=6 initial Ward partition.

This script does NOT select K=6 as final and does NOT run regime-restricted GWR.
It only inspects the K=6 partition already produced by the current exploratory
pipeline:

BasicGWR local slopes -> column-wise z-score -> Queen-constrained Ward -> K=6.

Outputs emphasize whether the six initial regions are spatially connected and
whether their local-relationship fingerprints are internally coherent.
"""

from __future__ import annotations

from pathlib import Path
import json

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[2]
COEF_CSV = ROOT / "results" / "real_data" / "georgia_gwr_baseline" / "basicgwr_local_coefficients.csv"
LABEL_CSV = ROOT / "results" / "spatial" / "georgia_initial_ward_regimes" / "initial_regime_labels_k2_k15.csv"
EDGE_CSV = ROOT / "results" / "spatial" / "georgia_queen_w" / "edge_list.csv"
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
OUT = ROOT / "results" / "spatial" / "georgia_initial_ward_regimes" / "k6_focus"

SLOPES = ["PctFB", "PctBlack", "PctRural"]


def _key(value) -> str:
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value)


def _standardize(values: np.ndarray) -> np.ndarray:
    mean = np.mean(values, axis=0)
    scale = np.std(values, axis=0, ddof=0)
    if np.any(scale <= 1e-12):
        raise RuntimeError(f"Near-constant slope encountered; std={scale.tolist()}")
    return (values - mean) / scale


def _build_neighbors(area_keys: list[str], edges: pd.DataFrame) -> tuple[dict[int, set[int]], list[tuple[int, int]]]:
    index = {key: i for i, key in enumerate(area_keys)}
    neighbors = {i: set() for i in range(len(area_keys))}
    pairs: set[tuple[int, int]] = set()

    for row in edges.itertuples(index=False):
        a = _key(row.AreaKey_i)
        b = _key(row.AreaKey_j)
        if a not in index or b not in index:
            raise RuntimeError(f"Unknown AreaKey in Queen edge list: {a}, {b}")
        i, j = index[a], index[b]
        neighbors[i].add(j)
        neighbors[j].add(i)
        pairs.add((min(i, j), max(i, j)))

    if len(pairs) != 431:
        raise RuntimeError(f"Expected 431 Queen edges; got {len(pairs)}")
    return neighbors, sorted(pairs)


def _is_connected(nodes: np.ndarray, neighbors: dict[int, set[int]]) -> bool:
    if nodes.size <= 1:
        return True
    allowed = set(int(v) for v in nodes)
    seen = {int(nodes[0])}
    stack = [int(nodes[0])]
    while stack:
        u = stack.pop()
        for v in neighbors[u]:
            if v in allowed and v not in seen:
                seen.add(v)
                stack.append(v)
    return len(seen) == len(allowed)


def _articulation_like_count(nodes: np.ndarray, neighbors: dict[int, set[int]]) -> int:
    """Count nodes whose removal disconnects the regime's induced Queen graph."""
    if nodes.size <= 2:
        return 0
    count = 0
    for node in nodes:
        remaining = nodes[nodes != node]
        if not _is_connected(remaining, neighbors):
            count += 1
    return count


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    coef = pd.read_csv(COEF_CSV)
    labels = pd.read_csv(LABEL_CSV)
    edges = pd.read_csv(EDGE_CSV)
    gdf = gpd.read_file(SHP)

    if len(coef) != 159 or len(labels) != 159 or len(gdf) != 159:
        raise RuntimeError("Expected 159 Georgia counties in all inputs")
    if "K6" not in labels.columns:
        raise RuntimeError("K6 labels are missing")

    coef["AreaKey_join"] = coef["AreaKey"].map(_key)
    labels["AreaKey_join"] = labels["AreaKey"].map(_key)
    table = coef.merge(labels[["AreaKey_join", "K6"]], on="AreaKey_join", how="left", validate="1:1")
    if table["K6"].isna().any():
        raise RuntimeError("Missing K6 labels after AreaKey join")

    table["K6"] = table["K6"].astype(int)
    if sorted(table["K6"].unique().tolist()) != [1, 2, 3, 4, 5, 6]:
        raise RuntimeError(f"Unexpected K6 IDs: {sorted(table['K6'].unique().tolist())}")

    features = _standardize(table[SLOPES].to_numpy(dtype=float))
    zcols = [f"z_{name}" for name in SLOPES]
    for j, col in enumerate(zcols):
        table[col] = features[:, j]

    area_keys = table["AreaKey_join"].tolist()
    neighbors, edge_pairs = _build_neighbors(area_keys, edges)

    summary_rows = []
    centers = []
    for regime in range(1, 7):
        nodes = np.flatnonzero(table["K6"].to_numpy() == regime)
        values = features[nodes]
        center = values.mean(axis=0)
        centers.append(center)
        distances = np.linalg.norm(values - center, axis=1)
        wcss = float(np.sum((values - center) ** 2))
        internal_edges = sum(i in set(nodes) and j in set(nodes) for i, j in edge_pairs)
        boundary_edges = sum((i in set(nodes)) != (j in set(nodes)) for i, j in edge_pairs)

        row = {
            "regime": regime,
            "n": int(nodes.size),
            "connected": bool(_is_connected(nodes, neighbors)),
            "articulation_like_nodes": int(_articulation_like_count(nodes, neighbors)),
            "internal_queen_edges": int(internal_edges),
            "incident_boundary_edges": int(boundary_edges),
            "wcss_contribution": wcss,
            "mean_distance_to_center": float(distances.mean()),
            "max_distance_to_center": float(distances.max()),
        }
        for j, slope in enumerate(SLOPES):
            row[f"mean_z_{slope}"] = float(center[j])
            row[f"std_z_{slope}"] = float(values[:, j].std(ddof=0))
            row[f"mean_raw_{slope}"] = float(table.iloc[nodes][slope].mean())
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "k6_regime_summary.csv", index=False)

    membership_cols = ["AreaKey", "K6", *SLOPES, *zcols]
    table[membership_cols].sort_values(["K6", "AreaKey"]).to_csv(
        OUT / "k6_regime_membership.csv", index=False
    )

    pair_rows = []
    for a in range(1, 7):
        for b in range(a + 1, 7):
            count = 0
            for i, j in edge_pairs:
                ri = int(table.iloc[i]["K6"])
                rj = int(table.iloc[j]["K6"])
                if {ri, rj} == {a, b}:
                    count += 1
            if count:
                pair_rows.append({"regime_a": a, "regime_b": b, "boundary_edges": count})
    pd.DataFrame(pair_rows).to_csv(OUT / "k6_regime_boundary_pairs.csv", index=False)

    shp_key = next((c for c in gdf.columns if str(c).lower() == "areakey"), None)
    if shp_key is None:
        raise RuntimeError("AreaKey not found in Georgia shapefile")
    gdf["AreaKey_join"] = gdf[shp_key].map(_key)
    mapped = gdf[["AreaKey_join", "geometry"]].merge(
        table[["AreaKey_join", "K6"]], on="AreaKey_join", how="left", validate="1:1"
    )
    if mapped["K6"].isna().any():
        raise RuntimeError("Missing K6 polygon labels")

    centers_arr = np.asarray(centers, dtype=float)
    heat_limit = float(np.max(np.abs(centers_arr)))
    heat_norm = TwoSlopeNorm(vmin=-heat_limit, vcenter=0.0, vmax=heat_limit)

    fig, (ax_map, ax_heat) = plt.subplots(
        1, 2, figsize=(15.5, 7.4), gridspec_kw={"width_ratios": [1.35, 1.0]}
    )

    mapped.plot(
        column="K6",
        categorical=True,
        cmap="tab10",
        ax=ax_map,
        edgecolor="black",
        linewidth=0.35,
        legend=True,
        legend_kwds={"title": "Initial regime"},
    )
    ax_map.set_title("Georgia exploratory K=6 initial partition")
    ax_map.set_aspect("equal")
    ax_map.set_axis_off()

    image = ax_heat.imshow(centers_arr, cmap="coolwarm", norm=heat_norm, aspect="auto")
    ax_heat.set_xticks(range(len(SLOPES)), SLOPES)
    ax_heat.set_yticks(range(6), [f"Regime {r}" for r in range(1, 7)])
    ax_heat.set_title("Mean standardized local-slope fingerprint")
    for i in range(6):
        for j in range(len(SLOPES)):
            ax_heat.text(j, i, f"{centers_arr[i, j]:.2f}", ha="center", va="center", fontsize=10)
    cbar = fig.colorbar(image, ax=ax_heat, shrink=0.78)
    cbar.set_label("Mean standardized BasicGWR slope")

    fig.suptitle(
        "K=6 is an exploratory working candidate, not a selected final regime count",
        fontsize=14,
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(OUT / "georgia_k6_initial_partition.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "georgia_k6_initial_partition.svg", bbox_inches="tight")
    plt.close(fig)

    global_boundary_edges = int(sum(table.iloc[i]["K6"] != table.iloc[j]["K6"] for i, j in edge_pairs))
    payload = {
        "status": "exploratory_working_candidate_only",
        "K": 6,
        "n_counties": 159,
        "regime_sizes": summary["n"].astype(int).tolist(),
        "all_regimes_connected": bool(summary["connected"].all()),
        "total_wcss": float(summary["wcss_contribution"].sum()),
        "queen_boundary_edges": global_boundary_edges,
        "queen_boundary_fraction": global_boundary_edges / len(edge_pairs),
        "fingerprint_terms": SLOPES,
        "intercept_used": False,
        "coordinates_used": False,
        "boundary_aware_refit_applied": False,
        "final_K_selected": False,
        "note": "Focused diagnostic only. K=6 remains an exploratory working candidate.",
    }
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(json.dumps(payload, indent=2))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
