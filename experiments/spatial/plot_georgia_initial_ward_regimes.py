"""Generate initial GR-GWR regime candidates from BasicGWR slope fingerprints.

This is an exploratory *initial segmentation* diagnostic, not the final GR-GWR
partition.  It intentionally follows the currently agreed clean construction:

1. use the repository BasicGWR local coefficient table already generated for
   the Georgia baseline;
2. exclude the intercept from the partition feature vector;
3. z-standardize the three local slopes (PctFB, PctBlack, PctRural);
4. use Queen contiguity only as the spatial connectivity constraint;
5. do NOT add coordinates to the feature vector;
6. fit spatially constrained Ward partitions for K=2,...,15.

The script persists all candidate labels and simple diagnostics so we can look
at the partitions before deciding how K should be selected or how later
boundary-aware refinement should work.
"""

from __future__ import annotations

from pathlib import Path
import json

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.cluster import AgglomerativeClustering

ROOT = Path(__file__).resolve().parents[2]
COEF_CSV = ROOT / "results" / "real_data" / "georgia_gwr_baseline" / "basicgwr_local_coefficients.csv"
EDGE_CSV = ROOT / "results" / "spatial" / "georgia_queen_w" / "edge_list.csv"
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
OUT = ROOT / "results" / "spatial" / "georgia_initial_ward_regimes"

SLOPES = ["PctFB", "PctBlack", "PctRural"]
K_VALUES = list(range(2, 16))


def _key(value) -> str:
    """Normalize AreaKey values so CSV/shapefile joins are stable."""
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


def _build_connectivity(area_keys: list[str], edges: pd.DataFrame) -> tuple[csr_matrix, list[tuple[int, int]]]:
    index = {key: i for i, key in enumerate(area_keys)}
    n = len(area_keys)
    graph = np.zeros((n, n), dtype=np.uint8)
    edge_pairs: list[tuple[int, int]] = []

    for row in edges.itertuples(index=False):
        a = _key(row.AreaKey_i)
        b = _key(row.AreaKey_j)
        if a not in index or b not in index:
            raise RuntimeError(f"Queen edge refers to unknown AreaKey: {a}, {b}")
        i, j = index[a], index[b]
        if i == j:
            raise RuntimeError("Self edge encountered in Queen edge list")
        graph[i, j] = 1
        graph[j, i] = 1
        edge_pairs.append((min(i, j), max(i, j)))

    edge_pairs = sorted(set(edge_pairs))
    if len(edge_pairs) != 431:
        raise RuntimeError(f"Expected 431 undirected Queen edges; got {len(edge_pairs)}")
    np.fill_diagonal(graph, 0)
    return csr_matrix(graph), edge_pairs


def _relabel_by_x_centroid(labels: np.ndarray, coords_x: np.ndarray) -> np.ndarray:
    """Give regimes deterministic left-to-right IDs for easier visual comparison."""
    old = np.unique(labels)
    order = sorted(old, key=lambda r: float(np.mean(coords_x[labels == r])))
    mapping = {int(r): new for new, r in enumerate(order)}
    return np.asarray([mapping[int(r)] for r in labels], dtype=int)


def _cluster_is_connected(nodes: np.ndarray, adjacency: csr_matrix) -> bool:
    if nodes.size <= 1:
        return True
    allowed = set(int(v) for v in nodes)
    seen = {int(nodes[0])}
    stack = [int(nodes[0])]
    while stack:
        u = stack.pop()
        neighbours = adjacency.indices[adjacency.indptr[u] : adjacency.indptr[u + 1]]
        for v in neighbours:
            v = int(v)
            if v in allowed and v not in seen:
                seen.add(v)
                stack.append(v)
    return len(seen) == len(allowed)


def _diagnostics(
    features: np.ndarray,
    labels: np.ndarray,
    edge_pairs: list[tuple[int, int]],
    adjacency: csr_matrix,
) -> dict[str, float | int | bool | str]:
    k = int(np.unique(labels).size)
    sizes = np.bincount(labels, minlength=k)

    wcss = 0.0
    connected = True
    for r in range(k):
        nodes = np.flatnonzero(labels == r)
        center = np.mean(features[nodes], axis=0)
        wcss += float(np.sum((features[nodes] - center) ** 2))
        connected = connected and _cluster_is_connected(nodes, adjacency)

    boundary_edges = int(sum(labels[i] != labels[j] for i, j in edge_pairs))
    return {
        "K": k,
        "wcss": wcss,
        "boundary_edges": boundary_edges,
        "boundary_fraction": boundary_edges / len(edge_pairs),
        "min_regime_size": int(sizes.min()),
        "median_regime_size": float(np.median(sizes)),
        "max_regime_size": int(sizes.max()),
        "all_regimes_connected": bool(connected),
        "regime_sizes": ";".join(str(int(v)) for v in sizes),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    coef = pd.read_csv(COEF_CSV)
    edges = pd.read_csv(EDGE_CSV)
    gdf = gpd.read_file(SHP)

    if len(coef) != 159 or len(gdf) != 159:
        raise RuntimeError(f"Expected 159 counties; coefficients={len(coef)}, polygons={len(gdf)}")

    coef["AreaKey_join"] = coef["AreaKey"].map(_key)
    area_keys = coef["AreaKey_join"].tolist()
    if len(set(area_keys)) != 159:
        raise RuntimeError("AreaKey is not unique in coefficient table")

    features = _standardize(coef[SLOPES].to_numpy(dtype=float))
    adjacency, edge_pairs = _build_connectivity(area_keys, edges)

    # Use polygon centroid X only to assign stable display IDs after clustering;
    # it never enters the Ward feature vector or merge criterion.
    shp_key = next((c for c in gdf.columns if str(c).lower() == "areakey"), None)
    if shp_key is None:
        raise RuntimeError("AreaKey not found in Georgia shapefile")
    gdf["AreaKey_join"] = gdf[shp_key].map(_key)
    x_by_key = dict(zip(gdf["AreaKey_join"], gdf.geometry.centroid.x))
    coords_x = np.asarray([x_by_key[k] for k in area_keys], dtype=float)

    labels_by_k: dict[int, np.ndarray] = {}
    diag_rows: list[dict[str, float | int | bool | str]] = []

    # Include K=1 diagnostics so the Ward merge cost K -> K-1 can be computed.
    labels_one = np.zeros(len(coef), dtype=int)
    diag_one = _diagnostics(features, labels_one, edge_pairs, adjacency)
    diagnostics_all = {1: diag_one}

    for k in K_VALUES:
        model = AgglomerativeClustering(
            n_clusters=k,
            linkage="ward",
            connectivity=adjacency,
            compute_full_tree=True,
        )
        labels = model.fit_predict(features)
        labels = _relabel_by_x_centroid(labels, coords_x)
        labels_by_k[k] = labels
        row = _diagnostics(features, labels, edge_pairs, adjacency)
        diagnostics_all[k] = row
        diag_rows.append(row)

    # Add the increase in WCSS when reducing K by one (the corresponding Ward merge cost).
    for row in diag_rows:
        k = int(row["K"])
        row["merge_cost_to_K_minus_1"] = float(diagnostics_all[k - 1]["wcss"] - row["wcss"])

    diag = pd.DataFrame(diag_rows)
    diag.to_csv(OUT / "ward_partition_diagnostics.csv", index=False)

    labels_table = pd.DataFrame({"AreaKey": coef["AreaKey"]})
    for k in K_VALUES:
        labels_table[f"K{k}"] = labels_by_k[k] + 1  # human-facing IDs start at 1
    labels_table.to_csv(OUT / "initial_regime_labels_k2_k15.csv", index=False)

    # Join all K labels to polygons once for plotting.
    mapped = gdf[["AreaKey_join", "geometry"]].copy()
    join = pd.DataFrame({"AreaKey_join": area_keys})
    for k in K_VALUES:
        join[f"K{k}"] = labels_by_k[k] + 1
    mapped = mapped.merge(join, on="AreaKey_join", how="left", validate="1:1")

    fig, axes = plt.subplots(4, 4, figsize=(19, 17))
    axes_flat = list(axes.flat)
    for ax, k in zip(axes_flat, K_VALUES):
        mapped.plot(
            column=f"K{k}",
            categorical=True,
            cmap="tab20",
            ax=ax,
            edgecolor="black",
            linewidth=0.25,
            legend=False,
        )
        row = diag.loc[diag["K"] == k].iloc[0]
        ax.set_title(
            f"K={k}  sizes {int(row['min_regime_size'])}-{int(row['max_regime_size'])}",
            fontsize=11,
        )
        ax.set_aspect("equal")
        ax.set_axis_off()

    for ax in axes_flat[len(K_VALUES):]:
        ax.set_axis_off()

    fig.suptitle(
        "Georgia initial regimes: Queen-constrained Ward on standardized BasicGWR slopes",
        fontsize=16,
        y=0.995,
    )
    fig.tight_layout()
    fig.savefig(OUT / "georgia_initial_ward_regimes_k2_k15.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "georgia_initial_ward_regimes_k2_k15.svg", bbox_inches="tight")
    plt.close(fig)

    # Diagnostics figure: no K is selected here; this only exposes the trade-off.
    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    ax.plot(diag["K"], diag["wcss"], marker="o")
    ax.set_xlabel("Number of regimes (K)")
    ax.set_ylabel("Within-regime coefficient SSE (WCSS)")
    ax.set_title("Initial Ward partition fit as K increases")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "ward_wcss_curve.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "ward_wcss_curve.svg", bbox_inches="tight")
    plt.close(fig)

    summary = {
        "n_counties": 159,
        "queen_edges": len(edge_pairs),
        "candidate_K": K_VALUES,
        "fingerprint_terms": SLOPES,
        "intercept_used": False,
        "coordinates_used_in_clustering_features": False,
        "spatial_constraint": "Queen contiguity",
        "clustering": "Ward agglomerative clustering constrained by Queen connectivity",
        "all_candidate_regimes_connected": bool(diag["all_regimes_connected"].all()),
        "note": "Exploratory initial partitions only; no final K selected and no boundary-aware refit applied.",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(diag.to_string(index=False))
    print(json.dumps(summary, indent=2))
    print(f"Generated: {(OUT / 'georgia_initial_ward_regimes_k2_k15.png').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
