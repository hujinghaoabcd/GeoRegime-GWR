"""Aggregate Simulation 3 scenario JSONs and draw phase-diagram summaries."""
from __future__ import annotations
import argparse
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def matrix(df, value):
    p = df.pivot(index="noise_sigma", columns="delta", values=value).sort_index(ascending=False)
    return p


def heatmap(df, value, title, output, center_zero=False):
    p = matrix(df, value)
    arr = p.to_numpy(dtype=float)
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    if center_zero:
        lim = np.nanmax(np.abs(arr))
        im = ax.imshow(arr, aspect="auto", cmap="coolwarm", vmin=-lim, vmax=lim)
    else:
        im = ax.imshow(arr, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(p.columns)), [f"{v:g}" for v in p.columns])
    ax.set_yticks(range(len(p.index)), [f"{v:g}" for v in p.index])
    ax.set_xlabel("Boundary strength δ")
    ax.set_ylabel("Noise σ")
    ax.set_title(title)
    for r in range(arr.shape[0]):
        for c in range(arr.shape[1]):
            v = arr[r,c]
            ax.text(c, r, "NA" if not np.isfinite(v) else f"{v:.3f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.82)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    inp = Path(args.input); out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    files = sorted(inp.rglob("scenario_*.json"))
    if len(files) != 20:
        raise RuntimeError(f"Expected 20 scenario files, found {len(files)}")
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in files]
    df = pd.DataFrame(rows).sort_values(["noise_sigma","delta"]).reset_index(drop=True)

    df["gr_minus_gwr_coef_rmse"] = df["grgwr_coef_rmse"] - df["gwr_coef_rmse"]
    df["gr_minus_mgwr_coef_rmse"] = df["grgwr_coef_rmse"] - df["mgwr_coef_rmse"]
    df["gr_vs_best_smoother_coef_rmse"] = df["grgwr_coef_rmse"] - df[["gwr_coef_rmse","mgwr_coef_rmse"]].min(axis=1)
    df["gr_oracle_gap_coef_rmse"] = df["grgwr_coef_rmse"] - df["oracle_coef_rmse"]
    df["gr_minus_gwr_fit_rmse"] = df["grgwr_fit_rmse"] - df["gwr_fit_rmse"]
    df["gr_minus_mgwr_fit_rmse"] = df["grgwr_fit_rmse"] - df["mgwr_fit_rmse"]
    df.to_csv(out / "phase_results.csv", index=False)

    heatmap(df, "gr_minus_gwr_coef_rmse", "GR-GWR minus GWR coefficient RMSE (negative is better)", out / "coef_rmse_gr_minus_gwr", True)
    heatmap(df, "gr_minus_mgwr_coef_rmse", "GR-GWR minus MGWR coefficient RMSE (negative is better)", out / "coef_rmse_gr_minus_mgwr", True)
    heatmap(df, "gr_vs_best_smoother_coef_rmse", "GR-GWR minus best of GWR/MGWR coefficient RMSE", out / "coef_rmse_gr_minus_best", True)
    heatmap(df[df.delta>0], "refined_boundary_f1", "Recovered boundary F1", out / "boundary_f1", False)
    heatmap(df[df.delta>0], "refined_ari", "Recovered regime ARI", out / "regime_ari", False)
    heatmap(df[df.delta>0], "grgwr_jump_ratio", "GR-GWR true boundary jump recovery ratio", out / "jump_recovery", False)

    positive = df[df.delta > 0].copy()
    win_gwr = int(np.sum(positive.grgwr_coef_rmse < positive.gwr_coef_rmse))
    win_mgwr = int(np.sum(positive.grgwr_coef_rmse < positive.mgwr_coef_rmse))
    win_best = int(np.sum(positive.grgwr_coef_rmse < positive[["gwr_coef_rmse","mgwr_coef_rmse"]].min(axis=1)))
    by_delta = positive.groupby("delta").agg(
        mean_boundary_f1=("refined_boundary_f1","mean"),
        mean_ari=("refined_ari","mean"),
        mean_gr_minus_gwr_coef_rmse=("gr_minus_gwr_coef_rmse","mean"),
        mean_gr_minus_mgwr_coef_rmse=("gr_minus_mgwr_coef_rmse","mean"),
        mean_gr_oracle_gap=("gr_oracle_gap_coef_rmse","mean"),
        mean_jump_recovery=("grgwr_jump_ratio","mean"),
    ).reset_index()
    by_noise = positive.groupby("noise_sigma").agg(
        mean_boundary_f1=("refined_boundary_f1","mean"),
        mean_ari=("refined_ari","mean"),
        mean_gr_minus_gwr_coef_rmse=("gr_minus_gwr_coef_rmse","mean"),
        mean_gr_minus_mgwr_coef_rmse=("gr_minus_mgwr_coef_rmse","mean"),
        mean_gr_oracle_gap=("gr_oracle_gap_coef_rmse","mean"),
    ).reset_index()
    by_delta.to_csv(out / "summary_by_delta.csv", index=False)
    by_noise.to_csv(out / "summary_by_noise.csv", index=False)

    summary = {
        "status":"single_seed_5x4_boundary_noise_phase_prototype",
        "n_scenarios":int(len(df)),
        "positive_boundary_scenarios":int(len(positive)),
        "delta_values":sorted(df.delta.unique().tolist()),
        "noise_values":sorted(df.noise_sigma.unique().tolist()),
        "wins_vs_gwr_coef_rmse":win_gwr,
        "wins_vs_mgwr_coef_rmse":win_mgwr,
        "wins_vs_best_gwr_mgwr_coef_rmse":win_best,
        "mean_positive_boundary_f1":float(positive.refined_boundary_f1.mean()),
        "mean_positive_ari":float(positive.refined_ari.mean()),
        "mean_gr_oracle_gap_coef_rmse":float(positive.gr_oracle_gap_coef_rmse.mean()),
        "warning":"Single seed per cell. K=3 is fixed. delta=0 is intentionally retained as a forced-K stress cell and is not a valid boundary-recovery target. Final publication evidence requires paired Monte Carlo repetitions and K-selection/no-boundary controls."
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("\nBy delta:\n", by_delta.to_string(index=False))
    print("\nBy noise:\n", by_noise.to_string(index=False))

if __name__ == "__main__": main()
