"""Aggregate paired Monte Carlo chunks for representative GR-GWR scenarios."""
from __future__ import annotations

import argparse
from pathlib import Path
import json
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def ci95(values: pd.Series) -> tuple[float, float, float, int]:
    x = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    n = x.size
    if n == 0:
        return float("nan"), float("nan"), float("nan"), 0
    mean = float(np.mean(x))
    if n == 1:
        return mean, mean, mean, 1
    se = float(np.std(x, ddof=1) / math.sqrt(n))
    half = 1.96 * se
    return mean, mean - half, mean + half, n


def scenario_summary(name: str, g: pd.DataFrame) -> dict:
    out: dict[str, float | int | str | None] = {
        "scenario": name,
        "delta": float(g["delta"].dropna().iloc[0]),
        "noise_sigma": float(g["noise_sigma"].dropna().iloc[0]),
        "requested_reps": int(g.shape[0]),
        "successful_reps": int(g["gwr_fit_rmse"].notna().sum()),
        "mgwr_successful_reps": int(g["mgwr_fit_rmse"].notna().sum()),
    }
    metrics = [
        "gwr_fit_rmse", "mgwr_fit_rmse", "grgwr_fit_rmse", "oracle_fit_rmse",
        "gwr_coef_rmse", "mgwr_coef_rmse", "grgwr_coef_rmse", "oracle_coef_rmse",
        "initial_ari", "refined_ari", "initial_boundary_f1", "refined_boundary_f1",
        "gwr_jump_ratio", "mgwr_jump_ratio", "grgwr_jump_ratio", "oracle_jump_ratio",
        "gr_minus_gwr_fit_rmse", "gr_minus_mgwr_fit_rmse",
        "gr_minus_gwr_coef_rmse", "gr_minus_mgwr_coef_rmse",
        "runtime_gwr_seconds", "runtime_mgwr_seconds", "runtime_gr_total_seconds",
        "runtime_gr_after_pilot_seconds", "refinement_iterations",
    ]
    for metric in metrics:
        if metric not in g.columns:
            continue
        mean, lo, hi, n = ci95(g[metric])
        out[f"{metric}_mean"] = mean
        out[f"{metric}_ci95_low"] = lo
        out[f"{metric}_ci95_high"] = hi
        out[f"{metric}_n"] = n

    win_cols = [
        "gr_fit_better_than_gwr", "gr_fit_better_than_mgwr",
        "gr_coef_better_than_gwr", "gr_coef_better_than_mgwr",
    ]
    for col in win_cols:
        if col in g.columns:
            vals = pd.to_numeric(g[col], errors="coerce").dropna()
            out[f"{col}_rate"] = float(vals.mean()) if len(vals) else float("nan")
            out[f"{col}_n"] = int(len(vals))

    if "refinement_stop" in g.columns:
        stop = g["refinement_stop"].dropna().value_counts().to_dict()
        out["refinement_stop_counts"] = json.dumps({str(k): int(v) for k, v in stop.items()}, sort_keys=True)
    return out


def paired_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    comparisons = [
        ("fit_rmse", "GWR", "gr_minus_gwr_fit_rmse"),
        ("fit_rmse", "MGWR", "gr_minus_mgwr_fit_rmse"),
        ("coef_rmse", "GWR", "gr_minus_gwr_coef_rmse"),
        ("coef_rmse", "MGWR", "gr_minus_mgwr_coef_rmse"),
    ]
    for scenario, g in df.groupby("scenario", sort=False):
        for metric, comparator, col in comparisons:
            if col not in g.columns:
                continue
            mean, lo, hi, n = ci95(g[col])
            rows.append({
                "scenario": scenario,
                "metric": metric,
                "comparison": f"GR-GWR minus {comparator}",
                "mean_paired_difference": mean,
                "ci95_low": lo,
                "ci95_high": hi,
                "n": n,
                "grgwr_better_if_negative": True,
                "ci_excludes_zero": bool(np.isfinite(lo) and np.isfinite(hi) and (hi < 0 or lo > 0)),
            })
    return pd.DataFrame(rows)


def plot_paired(summary: pd.DataFrame, out: Path) -> None:
    scenarios = summary["scenario"].tolist()
    x = np.arange(len(scenarios))

    fig, ax = plt.subplots(figsize=(10, 5.5))
    for label, col in [
        ("GR-GWR vs GWR", "gr_minus_gwr_fit_rmse"),
        ("GR-GWR vs MGWR", "gr_minus_mgwr_fit_rmse"),
    ]:
        means = summary[f"{col}_mean"].to_numpy(float)
        lo = summary[f"{col}_ci95_low"].to_numpy(float)
        hi = summary[f"{col}_ci95_high"].to_numpy(float)
        err = np.vstack([means - lo, hi - means])
        ax.errorbar(x, means, yerr=err, marker="o", capsize=4, label=label)
    ax.axhline(0.0, linewidth=1)
    ax.set_xticks(x, scenarios, rotation=25, ha="right")
    ax.set_ylabel("Paired RMSE difference (GR-GWR - comparator)")
    ax.set_title("Paired Monte Carlo fit-RMSE differences with 95% CI")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "paired_fit_rmse_differences.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / "paired_fit_rmse_differences.svg", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    for label, col in [
        ("GR-GWR vs GWR", "gr_minus_gwr_coef_rmse"),
        ("GR-GWR vs MGWR", "gr_minus_mgwr_coef_rmse"),
    ]:
        means = summary[f"{col}_mean"].to_numpy(float)
        lo = summary[f"{col}_ci95_low"].to_numpy(float)
        hi = summary[f"{col}_ci95_high"].to_numpy(float)
        err = np.vstack([means - lo, hi - means])
        ax.errorbar(x, means, yerr=err, marker="o", capsize=4, label=label)
    ax.axhline(0.0, linewidth=1)
    ax.set_xticks(x, scenarios, rotation=25, ha="right")
    ax.set_ylabel("Paired slope-coefficient RMSE difference")
    ax.set_title("Paired Monte Carlo coefficient-recovery differences with 95% CI")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "paired_coefficient_rmse_differences.png", dpi=300, bbox_inches="tight")
    fig.savefig(out / "paired_coefficient_rmse_differences.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    inp = Path(args.input)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    files = sorted(inp.rglob("mc_*_chunk*.csv"))
    if not files:
        raise RuntimeError(f"No Monte Carlo chunk CSVs found under {inp}")
    frames = []
    for f in files:
        df = pd.read_csv(f)
        # mc_<scenario>_chunkNN.csv; scenario may contain underscores.
        stem = f.stem
        scenario = stem[len("mc_"):].rsplit("_chunk", 1)[0]
        df["scenario"] = scenario
        df["source_file"] = f.name
        frames.append(df)
    all_df = pd.concat(frames, ignore_index=True)
    all_df = all_df.sort_values(["scenario", "seed"]).reset_index(drop=True)
    all_df.to_csv(out / "all_replicates.csv", index=False)

    summary_rows = [scenario_summary(name, g) for name, g in all_df.groupby("scenario", sort=False)]
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out / "scenario_summary.csv", index=False)

    paired = paired_table(all_df)
    paired.to_csv(out / "paired_differences.csv", index=False)
    plot_paired(summary, out)

    payload = {
        "status": "paired_monte_carlo_representative_scenarios",
        "total_rows": int(len(all_df)),
        "scenarios": int(all_df["scenario"].nunique()),
        "scenario_names": summary["scenario"].tolist(),
        "target_reps_per_scenario": 100,
        "paired_design": "within each seed, all compared models use the same generated X/y realization",
        "K_fixed": 3,
        "true_boundaries_hidden_from_estimated_grgwr": True,
        "lambda_working_baseline": 0.5,
        "primary_interpretation": "Coefficient recovery and boundary recovery are primary; in-sample fit is secondary. The delta=0 forced-K scenario is a stress control, not a valid boundary-recovery target.",
        "scenario_summary": summary.to_dict(orient="records"),
    }
    (out / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    print("\nPaired differences:\n", paired.to_string(index=False))


if __name__ == "__main__":
    main()
