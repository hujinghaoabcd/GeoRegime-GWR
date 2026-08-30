"""Visualize how K=6 labels change after the current unified GR-GWR refinement.

Outputs:
- changed_counties.csv: the 35 counties whose final regime differs from the initial Ward regime;
- changed_transition_counts.csv: initial->final move counts;
- regime_size_change.csv: initial/final regime sizes;
- georgia_k6_refinement_changes.png/.svg: initial map, final map, and transition matrix.

This is exploratory visualization only; it does not validate the final partition.
"""

from __future__ import annotations

from pathlib import Path
import json

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[2]
SHP = ROOT / "data" / "raw" / "georgia" / "G_utm.shp"
RESULTS = ROOT / "results" / "real_data" / "georgia_k6_unified_label_refinement"
COUNTY_CSV = RESULTS / "final_county_results.csv"


def _key(value) -> str:
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value)


def main() -> None:
    county = pd.read_csv(COUNTY_CSV)
    gdf = gpd.read_file(SHP)
    if len(county) != 159 or len(gdf) != 159:
        raise RuntimeError("Expected 159 Georgia counties")

    county["AreaKey_join"] = county["AreaKey"].map(_key)
    shp_key = next((c for c in gdf.columns if str(c).lower() == "areakey"), None)
    if shp_key is None:
        raise RuntimeError("AreaKey not found in Georgia shapefile")
    gdf["AreaKey_join"] = gdf[shp_key].map(_key)

    mapped = gdf[["AreaKey_join", "geometry"]].merge(
        county[["AreaKey_join", "AreaKey", "initial_regime", "final_regime", "regime_changed"]],
        on="AreaKey_join",
        how="left",
        validate="1:1",
    )
    if mapped[["initial_regime", "final_regime"]].isna().any().any():
        raise RuntimeError("Missing regime labels after shapefile join")

    mapped["initial_regime"] = mapped["initial_regime"].astype(int)
    mapped["final_regime"] = mapped["final_regime"].astype(int)
    mapped["regime_changed"] = mapped["regime_changed"].astype(bool)

    changed = mapped.loc[mapped["regime_changed"], ["AreaKey", "initial_regime", "final_regime"]].copy()
    changed["transition"] = changed["initial_regime"].astype(str) + "→" + changed["final_regime"].astype(str)
    changed = changed.sort_values(["initial_regime", "final_regime", "AreaKey"])
    changed.to_csv(RESULTS / "changed_counties.csv", index=False)

    matrix = pd.crosstab(changed["initial_regime"], changed["final_regime"]).reindex(
        index=range(1, 7), columns=range(1, 7), fill_value=0
    )
    matrix.index.name = "from_regime"
    matrix.columns = [f"to_{c}" for c in matrix.columns]
    matrix.to_csv(RESULTS / "changed_transition_counts.csv")

    size_rows = []
    for r in range(1, 7):
        n0 = int((mapped["initial_regime"] == r).sum())
        n1 = int((mapped["final_regime"] == r).sum())
        size_rows.append({"regime": r, "initial_n": n0, "final_n": n1, "delta_n": n1 - n0})
    size_df = pd.DataFrame(size_rows)
    size_df.to_csv(RESULTS / "regime_size_change.csv", index=False)

    transition_long = (
        changed.groupby(["initial_regime", "final_regime"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
        .sort_values(["count", "initial_regime", "final_regime"], ascending=[False, True, True])
    )
    transition_long.to_csv(RESULTS / "changed_transition_summary.csv", index=False)

    base = plt.get_cmap("tab10")
    cmap = ListedColormap([base(i) for i in range(6)])
    norm = BoundaryNorm(np.arange(0.5, 7.5, 1.0), cmap.N)

    fig, axes = plt.subplots(1, 3, figsize=(19, 7.4), gridspec_kw={"width_ratios": [1.12, 1.12, 0.9]})
    ax0, ax1, ax2 = axes

    mapped.plot(column="initial_regime", cmap=cmap, norm=norm, ax=ax0, edgecolor="black", linewidth=0.3)
    ax0.set_title("Initial K=6 Ward partition\n53, 19, 17, 22, 27, 21")
    ax0.set_axis_off()
    ax0.set_aspect("equal")

    mapped.plot(column="final_regime", cmap=cmap, norm=norm, ax=ax1, edgecolor="black", linewidth=0.3)
    mapped.loc[mapped["regime_changed"]].boundary.plot(ax=ax1, color="black", linewidth=1.25)
    final_sizes = size_df["final_n"].tolist()
    ax1.set_title("After 8 accepted refinement iterations\n" + ", ".join(map(str, final_sizes)) + " (thick outlines = changed)")
    ax1.set_axis_off()
    ax1.set_aspect("equal")

    matrix_values = pd.crosstab(changed["initial_regime"], changed["final_regime"]).reindex(
        index=range(1, 7), columns=range(1, 7), fill_value=0
    ).to_numpy()
    im = ax2.imshow(matrix_values, cmap="Blues", aspect="equal")
    ax2.set_xticks(range(6), [f"R{r}" for r in range(1, 7)])
    ax2.set_yticks(range(6), [f"R{r}" for r in range(1, 7)])
    ax2.set_xlabel("Final regime")
    ax2.set_ylabel("Initial regime")
    ax2.set_title("35 changed counties: transition counts")
    for i in range(6):
        for j in range(6):
            value = int(matrix_values[i, j])
            if value:
                ax2.text(j, i, str(value), ha="center", va="center", fontsize=11)
    fig.colorbar(im, ax=ax2, shrink=0.72, label="County count")

    legend = [Patch(facecolor=cmap(r - 1), edgecolor="black", label=f"Regime {r}") for r in range(1, 7)]
    fig.legend(handles=legend, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.42, 0.01))
    fig.suptitle("Georgia K=6 refinement: where labels changed and the final spatial partition", fontsize=15, y=0.99)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    fig.savefig(RESULTS / "georgia_k6_refinement_changes.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULTS / "georgia_k6_refinement_changes.svg", bbox_inches="tight")
    plt.close(fig)

    payload = {
        "changed_counties": int(len(changed)),
        "unchanged_counties": int(159 - len(changed)),
        "initial_sizes": size_df["initial_n"].astype(int).tolist(),
        "final_sizes": size_df["final_n"].astype(int).tolist(),
        "size_changes": size_df["delta_n"].astype(int).tolist(),
        "transition_counts": [
            {"from": int(r.initial_regime), "to": int(r.final_regime), "count": int(r.count)}
            for r in transition_long.itertuples(index=False)
        ],
        "note": "Exploratory same-sample refinement visualization; not final partition validation.",
    }
    (RESULTS / "label_change_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print("\nChanged counties:\n", changed.to_string(index=False))


if __name__ == "__main__":
    main()
