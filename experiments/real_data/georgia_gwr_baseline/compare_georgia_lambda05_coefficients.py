"""Compare local coefficients for GWR, MGWR, and the lambda=0.5 GR-GWR working baseline.

This is a diagnostic, not a final inference experiment.  The goal is to see
how regime-aware refinement changes local coefficients and whether those
changes/improvements concentrate near final regime boundaries.
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from georegime_gwr import BasicGWR, RegimeAwareGWR

DATA = ROOT / "data" / "raw" / "georgia" / "GData_utm.csv"
LABELS = ROOT / "results" / "real_data" / "georgia_k6_lambda_sensitivity" / "lambda_final_labels.csv"
MGWR = ROOT / "results" / "real_data" / "georgia_mgwr_benchmark" / "mgwr_county_results.csv"
EDGES = ROOT / "results" / "spatial" / "georgia_queen_w" / "edge_list.csv"
OUT = ROOT / "results" / "real_data" / "georgia_lambda05_coefficient_comparison"

PREDICTORS = ["PctFB", "PctBlack", "PctRural"]
TERMS = ["Intercept", *PREDICTORS]
LAMBDA_COLUMN = "lambda_0.5"


def _key(value) -> str:
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value)


def _metrics(y: np.ndarray, fitted: np.ndarray) -> dict[str, float]:
    residuals = y - fitted
    rss = float(residuals @ residuals)
    return {
        "rss": rss,
        "rmse": float(np.sqrt(np.mean(residuals**2))),
        "mae": float(np.mean(np.abs(residuals))),
        "r2": float(1.0 - rss / np.sum((y - np.mean(y)) ** 2)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(DATA)
    labels_df = pd.read_csv(LABELS)
    mgwr_df = pd.read_csv(MGWR)
    edge_df = pd.read_csv(EDGES)

    for frame in (df, labels_df, mgwr_df):
        frame["AreaKey_join"] = frame["AreaKey"].map(_key)

    df = df.merge(
        labels_df[["AreaKey_join", "initial", LAMBDA_COLUMN]],
        on="AreaKey_join", how="left", validate="1:1"
    )
    if len(df) != 159 or df[LAMBDA_COLUMN].isna().any():
        raise RuntimeError("Expected complete lambda=0.5 labels for 159 counties")

    mgwr_cols = ["AreaKey_join", *TERMS, "fitted_z", "residual_z"]
    mgwr_df = mgwr_df[mgwr_cols].copy()
    mgwr_df = mgwr_df.rename(columns={
        **{term: f"mgwr_{term}" for term in TERMS},
        "fitted_z": "mgwr_fitted_z",
        "residual_z": "mgwr_residual_z",
    })
    df = df.merge(mgwr_df, on="AreaKey_join", how="left", validate="1:1")

    X = df[PREDICTORS].to_numpy(dtype=float)
    y = df["PctBach"].to_numpy(dtype=float)
    coords = df[["X", "Y"]].to_numpy(dtype=float)
    Xz = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    yz = (y - y.mean()) / y.std(ddof=0)
    labels = df[LAMBDA_COLUMN].to_numpy(dtype=int)
    initial = df["initial"].to_numpy(dtype=int)

    gwr = BasicGWR(bandwidth="auto", kernel="bisquare", fit_intercept=True).fit(Xz, yz, coords)
    if int(gwr.bandwidth_) != 116:
        raise RuntimeError(f"Expected ordinary GWR k=116, got {gwr.bandwidth_}")

    gr = RegimeAwareGWR(
        bandwidth="regime_size", kernel="bisquare", fit_intercept=True
    ).fit(Xz, yz, coords, labels)

    n = len(df)
    adjacency = [set() for _ in range(n)]
    for row in edge_df.itertuples(index=False):
        i, j = int(row.i), int(row.j)
        adjacency[i].add(j)
        adjacency[j].add(i)
    boundary = np.array([
        any(labels[j] != labels[i] for j in adjacency[i]) for i in range(n)
    ], dtype=bool)

    county = pd.DataFrame({
        "AreaKey": df["AreaKey"],
        "initial_regime": initial,
        "final_regime_lambda05": labels,
        "changed_from_initial": initial != labels,
        "final_boundary_county": boundary,
        "gwr_fitted_z": gwr.fitted_values_,
        "gwr_residual_z": gwr.residuals_,
        "mgwr_fitted_z": df["mgwr_fitted_z"],
        "mgwr_residual_z": df["mgwr_residual_z"],
        "grgwr_fitted_z": gr.fitted_values_,
        "grgwr_residual_z": gr.residuals_,
    })
    for j, term in enumerate(TERMS):
        county[f"gwr_{term}"] = gwr.parameters_[:, j]
        county[f"mgwr_{term}"] = df[f"mgwr_{term}"].to_numpy(dtype=float)
        county[f"grgwr_{term}"] = gr.parameters_[:, j]

    county["abs_resid_improvement_gr_vs_gwr"] = (
        np.abs(county["gwr_residual_z"]) - np.abs(county["grgwr_residual_z"])
    )
    county["abs_resid_improvement_gr_vs_mgwr"] = (
        np.abs(county["mgwr_residual_z"]) - np.abs(county["grgwr_residual_z"])
    )
    county.to_csv(OUT / "county_coefficient_comparison.csv", index=False)

    pair_rows = []
    for term in TERMS:
        arrays = {
            "gwr": county[f"gwr_{term}"].to_numpy(dtype=float),
            "mgwr": county[f"mgwr_{term}"].to_numpy(dtype=float),
            "grgwr": county[f"grgwr_{term}"].to_numpy(dtype=float),
        }
        for a, b in [("gwr", "mgwr"), ("gwr", "grgwr"), ("mgwr", "grgwr")]:
            diff = arrays[b] - arrays[a]
            pair_rows.append({
                "term": term,
                "model_a": a,
                "model_b": b,
                "pearson_r": float(np.corrcoef(arrays[a], arrays[b])[0, 1]),
                "mean_abs_difference": float(np.mean(np.abs(diff))),
                "rmse_difference": float(np.sqrt(np.mean(diff**2))),
                "max_abs_difference": float(np.max(np.abs(diff))),
            })
    pair_df = pd.DataFrame(pair_rows)
    pair_df.to_csv(OUT / "coefficient_pair_summary.csv", index=False)

    boundary_rows = []
    for group_name, mask in [("boundary", boundary), ("interior", ~boundary)]:
        row = {
            "group": group_name,
            "n": int(np.sum(mask)),
            "gwr_mean_abs_residual": float(np.mean(np.abs(gwr.residuals_[mask]))),
            "mgwr_mean_abs_residual": float(np.mean(np.abs(df.loc[mask, "mgwr_residual_z"].to_numpy(dtype=float)))),
            "grgwr_mean_abs_residual": float(np.mean(np.abs(gr.residuals_[mask]))),
            "mean_abs_resid_improvement_gr_vs_gwr": float(np.mean(county.loc[mask, "abs_resid_improvement_gr_vs_gwr"])),
            "mean_abs_resid_improvement_gr_vs_mgwr": float(np.mean(county.loc[mask, "abs_resid_improvement_gr_vs_mgwr"])),
        }
        boundary_rows.append(row)
    boundary_df = pd.DataFrame(boundary_rows)
    boundary_df.to_csv(OUT / "boundary_residual_summary.csv", index=False)

    coef_boundary_rows = []
    for term in TERMS:
        for group_name, mask in [("boundary", boundary), ("interior", ~boundary)]:
            gwr_vals = county.loc[mask, f"gwr_{term}"].to_numpy(dtype=float)
            mgwr_vals = county.loc[mask, f"mgwr_{term}"].to_numpy(dtype=float)
            gr_vals = county.loc[mask, f"grgwr_{term}"].to_numpy(dtype=float)
            coef_boundary_rows.append({
                "term": term,
                "group": group_name,
                "n": int(np.sum(mask)),
                "mean_abs_gr_minus_gwr": float(np.mean(np.abs(gr_vals - gwr_vals))),
                "mean_abs_gr_minus_mgwr": float(np.mean(np.abs(gr_vals - mgwr_vals))),
                "mean_abs_mgwr_minus_gwr": float(np.mean(np.abs(mgwr_vals - gwr_vals))),
            })
    coef_boundary_df = pd.DataFrame(coef_boundary_rows)
    coef_boundary_df.to_csv(OUT / "coefficient_boundary_difference.csv", index=False)

    changed_rows = []
    changed_mask = initial != labels
    for group_name, mask in [("changed", changed_mask), ("unchanged", ~changed_mask)]:
        changed_rows.append({
            "group": group_name,
            "n": int(np.sum(mask)),
            "gwr_mean_abs_residual": float(np.mean(np.abs(gwr.residuals_[mask]))),
            "grgwr_mean_abs_residual": float(np.mean(np.abs(gr.residuals_[mask]))),
            "mean_abs_resid_improvement_gr_vs_gwr": float(np.mean(county.loc[mask, "abs_resid_improvement_gr_vs_gwr"])),
        })
    changed_df = pd.DataFrame(changed_rows)
    changed_df.to_csv(OUT / "changed_county_residual_summary.csv", index=False)

    summary = {
        "status": "exploratory_lambda05_coefficient_diagnostic",
        "lambda_working_baseline": 0.5,
        "n": 159,
        "final_regime_sizes": {str(r): int(np.sum(labels == r)) for r in np.unique(labels)},
        "changed_from_initial": int(np.sum(changed_mask)),
        "boundary_counties": int(np.sum(boundary)),
        "interior_counties": int(np.sum(~boundary)),
        "ordinary_gwr": _metrics(yz, gwr.fitted_values_),
        "mgwr": _metrics(yz, df["mgwr_fitted_z"].to_numpy(dtype=float)),
        "grgwr_lambda05": _metrics(yz, gr.fitted_values_),
        "boundary_residual_summary": boundary_df.to_dict(orient="records"),
        "changed_county_residual_summary": changed_df.to_dict(orient="records"),
        "interpretation_guard": (
            "All results are in-sample diagnostics. Final boundary labels are data-adaptive; "
            "boundary/interior contrasts are descriptive and not causal or out-of-sample evidence."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print("\nCoefficient pair summary:\n", pair_df.to_string(index=False))
    print("\nBoundary residual summary:\n", boundary_df.to_string(index=False))
    print("\nCoefficient boundary differences:\n", coef_boundary_df.to_string(index=False))


if __name__ == "__main__":
    main()
