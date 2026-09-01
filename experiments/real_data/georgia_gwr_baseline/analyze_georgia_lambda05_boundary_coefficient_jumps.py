"""Map and quantify coefficient discontinuities at final lambda=0.5 GR-GWR boundaries.

Exploratory diagnostic only.  This script does not re-select lambda or regimes.
It asks a narrower descriptive question: relative to ordinary GWR/MGWR, does the
working GR-GWR baseline show larger coefficient contrasts across final regime
boundaries than within regimes, as expected if cross-boundary smoothing is
being reduced?
"""

from __future__ import annotations

from pathlib import Path
import json

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
COUNTY = ROOT / "results" / "real_data" / "georgia_lambda05_coefficient_comparison" / "county_coefficient_comparison.csv"
EDGES = ROOT / "results" / "spatial" / "georgia_queen_w" / "edge_list.csv"
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
OUT = ROOT / "results" / "real_data" / "georgia_lambda05_boundary_coefficient_jumps"

TERMS = ["PctFB", "PctBlack", "PctRural"]
MODELS = ["gwr", "mgwr", "grgwr"]


def _key(value) -> str:
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value)


def _norm(values: np.ndarray):
    values = np.asarray(values, dtype=float)
    vmin = float(np.nanmin(values))
    vmax = float(np.nanmax(values))
    if vmin < 0.0 < vmax:
        return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    if np.isclose(vmin, vmax):
        eps = max(abs(vmin) * 1e-6, 1e-9)
        return Normalize(vmin=vmin - eps, vmax=vmax + eps)
    return Normalize(vmin=vmin, vmax=vmax)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    county = pd.read_csv(COUNTY)
    edges = pd.read_csv(EDGES)
    gdf = gpd.read_file(SHP)

    if len(county) != 159:
        raise RuntimeError(f"Expected 159 counties, got {len(county)}")
    if not {"i", "j"}.issubset(edges.columns):
        raise RuntimeError("Queen edge list must contain i and j columns")

    labels = county["final_regime_lambda05"].to_numpy(dtype=int)
    if sorted(np.unique(labels).tolist()) != [1, 2, 3, 4, 5, 6]:
        raise RuntimeError(f"Unexpected final regimes: {np.unique(labels).tolist()}")

    edge_rows: list[dict[str, float | int | str | bool]] = []
    for edge_id, row in enumerate(edges.itertuples(index=False)):
        i, j = int(row.i), int(row.j)
        boundary = bool(labels[i] != labels[j])
        base = {
            "edge_id": edge_id,
            "i": i,
            "j": j,
            "AreaKey_i": county.iloc[i]["AreaKey"],
            "AreaKey_j": county.iloc[j]["AreaKey"],
            "regime_i": int(labels[i]),
            "regime_j": int(labels[j]),
            "boundary_edge": boundary,
        }
        for term in TERMS:
            for model in MODELS:
                a = float(county.iloc[i][f"{model}_{term}"])
                b = float(county.iloc[j][f"{model}_{term}"])
                base[f"{model}_{term}_abs_jump"] = abs(a - b)
        edge_rows.append(base)

    edge_df = pd.DataFrame(edge_rows)
    edge_df.to_csv(OUT / "queen_edge_coefficient_jumps.csv", index=False)

    summary_rows = []
    for term in TERMS:
        for model in MODELS:
            col = f"{model}_{term}_abs_jump"
            boundary_vals = edge_df.loc[edge_df["boundary_edge"], col].to_numpy(dtype=float)
            interior_vals = edge_df.loc[~edge_df["boundary_edge"], col].to_numpy(dtype=float)
            summary_rows.append({
                "term": term,
                "model": model,
                "boundary_edges": int(boundary_vals.size),
                "interior_edges": int(interior_vals.size),
                "boundary_mean_abs_jump": float(np.mean(boundary_vals)),
                "boundary_median_abs_jump": float(np.median(boundary_vals)),
                "interior_mean_abs_jump": float(np.mean(interior_vals)),
                "interior_median_abs_jump": float(np.median(interior_vals)),
                "boundary_to_interior_mean_ratio": float(np.mean(boundary_vals) / np.mean(interior_vals)) if np.mean(interior_vals) > 0 else float("nan"),
            })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT / "edge_jump_summary.csv", index=False)

    amp_rows = []
    boundary_mask = edge_df["boundary_edge"].to_numpy(dtype=bool)
    for term in TERMS:
        gwr = edge_df.loc[boundary_mask, f"gwr_{term}_abs_jump"].to_numpy(dtype=float)
        mgwr = edge_df.loc[boundary_mask, f"mgwr_{term}_abs_jump"].to_numpy(dtype=float)
        gr = edge_df.loc[boundary_mask, f"grgwr_{term}_abs_jump"].to_numpy(dtype=float)
        amp_rows.append({
            "term": term,
            "boundary_edges": int(gr.size),
            "mean_gr_minus_gwr_jump": float(np.mean(gr - gwr)),
            "median_gr_minus_gwr_jump": float(np.median(gr - gwr)),
            "fraction_boundary_edges_gr_jump_gt_gwr": float(np.mean(gr > gwr)),
            "mean_gr_to_gwr_jump_ratio": float(np.mean(gr) / np.mean(gwr)) if np.mean(gwr) > 0 else float("nan"),
            "mean_gr_minus_mgwr_jump": float(np.mean(gr - mgwr)),
            "fraction_boundary_edges_gr_jump_gt_mgwr": float(np.mean(gr > mgwr)),
            "mean_gr_to_mgwr_jump_ratio": float(np.mean(gr) / np.mean(mgwr)) if np.mean(mgwr) > 0 else float("nan"),
        })
    amp_df = pd.DataFrame(amp_rows)
    amp_df.to_csv(OUT / "boundary_jump_amplification.csv", index=False)

    # Join county coefficients to polygons.
    shp_key = next((c for c in gdf.columns if str(c).lower() == "areakey"), None)
    if shp_key is None:
        raise RuntimeError("AreaKey not found in Georgia shapefile")
    gdf["AreaKey_join"] = gdf[shp_key].map(_key)
    county["AreaKey_join"] = county["AreaKey"].map(_key)
    map_df = gdf[["AreaKey_join", "geometry"]].merge(county, on="AreaKey_join", how="left", validate="1:1")
    if map_df["final_regime_lambda05"].isna().any():
        raise RuntimeError("Failed to join coefficient results to polygons")

    regime_outline = map_df[["final_regime_lambda05", "geometry"]].dissolve(by="final_regime_lambda05")

    # One 3-panel figure per slope: GWR, GR-GWR, and GR-GWR minus GWR.
    for term in TERMS:
        gwr_col = f"gwr_{term}"
        gr_col = f"grgwr_{term}"
        delta_col = f"delta_gr_minus_gwr_{term}"
        map_df[delta_col] = map_df[gr_col] - map_df[gwr_col]

        combined = np.concatenate([
            map_df[gwr_col].to_numpy(dtype=float),
            map_df[gr_col].to_numpy(dtype=float),
        ])
        shared_norm = _norm(combined)
        delta_norm = _norm(map_df[delta_col].to_numpy(dtype=float))

        fig, axes = plt.subplots(1, 3, figsize=(18, 6.2))
        for ax, col, title, norm in [
            (axes[0], gwr_col, f"GWR: {term}", shared_norm),
            (axes[1], gr_col, f"GR-GWR (lambda=0.5): {term}", shared_norm),
            (axes[2], delta_col, f"GR-GWR - GWR: {term}", delta_norm),
        ]:
            map_df.plot(
                column=col,
                cmap="coolwarm",
                norm=norm,
                edgecolor="white",
                linewidth=0.25,
                ax=ax,
                legend=True,
                legend_kwds={"shrink": 0.72},
            )
            regime_outline.boundary.plot(ax=ax, linewidth=1.15, edgecolor="black")
            ax.set_axis_off()
            ax.set_title(title)
        fig.suptitle(f"Georgia local coefficient comparison with final regime boundaries: {term}", fontsize=15)
        fig.tight_layout()
        fig.savefig(OUT / f"{term}_gwr_grgwr_delta_map.png", dpi=300, bbox_inches="tight")
        fig.savefig(OUT / f"{term}_gwr_grgwr_delta_map.svg", bbox_inches="tight")
        plt.close(fig)

    # Compact cross-boundary jump comparison chart.
    chart = amp_df.set_index("term")[["mean_gr_to_gwr_jump_ratio", "mean_gr_to_mgwr_jump_ratio"]]
    ax = chart.plot(kind="bar", figsize=(8.6, 5.4))
    ax.axhline(1.0, linewidth=1.0, linestyle="--")
    ax.set_ylabel("Mean cross-boundary jump ratio")
    ax.set_xlabel("")
    ax.set_title("Cross-boundary coefficient jump: GR-GWR relative to GWR/MGWR")
    ax.legend(["GR-GWR / GWR", "GR-GWR / MGWR"])
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(OUT / "cross_boundary_jump_ratio.png", dpi=300, bbox_inches="tight")
    plt.savefig(OUT / "cross_boundary_jump_ratio.svg", bbox_inches="tight")
    plt.close()

    summary = {
        "status": "exploratory_descriptive_boundary_coefficient_jump_diagnostic",
        "lambda_working_baseline": 0.5,
        "n_counties": 159,
        "queen_edges": int(len(edge_df)),
        "final_boundary_edges": int(edge_df["boundary_edge"].sum()),
        "within_regime_edges": int((~edge_df["boundary_edge"]).sum()),
        "boundary_jump_amplification": amp_df.to_dict(orient="records"),
        "edge_jump_summary": summary_df.to_dict(orient="records"),
        "interpretation_guard": (
            "The final regimes are data-adaptive and were learned from the same response data. "
            "These edge contrasts are descriptive diagnostics, not independent evidence that the inferred boundaries are true."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("\nBoundary jump amplification:\n", amp_df.to_string(index=False))
    print("\nEdge jump summary:\n", summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
