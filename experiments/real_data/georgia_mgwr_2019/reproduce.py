"""Reproduce the Georgia GWR/MGWR example from Oshan et al. (2019).

This script intentionally uses the external ``mgwr`` package rather than the
research implementation in this repository.  Its purpose is to establish a
trusted baseline before GR-GWR experiments are introduced.

Reference model specification
-----------------------------
Dependent variable:
    PctBach
Explanatory variables:
    PctFB, PctBlack, PctRural
Coordinates:
    projected X, Y
Pre-processing:
    z-score y and each X column using population mean/std (ddof=0)
Kernel:
    adaptive bisquare
Bandwidth criterion:
    AICc

The historical PySAL notebook reports:
    GWR bandwidth = 117
    MGWR bandwidths = [92, 101, 136, 158]

Run from repository root:
    python experiments/real_data/georgia_mgwr_2019/reproduce.py --vendor-data
"""

from __future__ import annotations

import argparse
import json
import shutil
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd
import libpysal as ps
from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW


ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = ROOT / "data" / "raw" / "georgia"
RESULT_DIR = ROOT / "results" / "reproduction" / "georgia_mgwr_2019"

VARIABLE_NAMES = ["Intercept", "PctFB", "PctBlack", "PctRural"]

REFERENCE = {
    "source": "PySAL MGWR_Georgia_example / Oshan et al. (2019)",
    "n": 159,
    "gwr": {
        "bandwidth": 117.0,
        "rss": 51.186,
        "enp": 11.805,
        "aic": 296.616,
        "aicc": 299.051,
        "bic": 335.913,
        "r2": 0.678,
        "adj_r2": 0.652,
    },
    "mgwr": {
        "bandwidths": [92.0, 101.0, 136.0, 158.0],
        "rss": 50.899,
        "enp": 11.368,
        "aic": 294.849,
        "aicc": 297.120,
        "bic": 332.806,
        "r2": 0.680,
        "adj_r2": 0.655,
    },
}


def _scalar(value) -> float:
    arr = np.asarray(value).reshape(-1)
    return float(arr[0])


def vendor_georgia_data() -> Path:
    """Copy the canonical libpysal Georgia example into this research repo."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    csv_src = Path(ps.examples.get_path("GData_utm.csv"))
    csv_dst = DATA_DIR / "GData_utm.csv"
    shutil.copy2(csv_src, csv_dst)

    # Preserve county polygons for later spatial-weight and mapping experiments.
    shp_src = Path(ps.examples.get_path("G_utm.shp"))
    for src in shp_src.parent.glob("G_utm.*"):
        if src.is_file():
            shutil.copy2(src, DATA_DIR / src.name)

    return csv_dst


def locate_data(vendor: bool) -> Path:
    local = DATA_DIR / "GData_utm.csv"
    if vendor:
        return vendor_georgia_data()
    if local.exists():
        return local
    return Path(ps.examples.get_path("GData_utm.csv"))


def parameter_summary(params: np.ndarray) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for j, name in enumerate(VARIABLE_NAMES):
        values = params[:, j]
        out[name] = {
            "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=0)),
            "min": float(np.min(values)),
            "median": float(np.median(values)),
            "max": float(np.max(values)),
        }
    return out


def fit_models(data_path: Path):
    data = pd.read_csv(data_path)

    y = data["PctBach"].to_numpy(dtype=float).reshape((-1, 1))
    X = data[["PctFB", "PctBlack", "PctRural"]].to_numpy(dtype=float)
    coords = list(zip(data["X"].to_numpy(dtype=float), data["Y"].to_numpy(dtype=float)))

    # Exactly match the published notebook preprocessing (ddof=0).
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    y = (y - y.mean(axis=0)) / y.std(axis=0)

    gwr_selector = Sel_BW(coords, y, X)
    gwr_bw = gwr_selector.search(bw_min=2)
    gwr_results = GWR(coords, y, X, gwr_bw).fit()

    mgwr_selector = Sel_BW(coords, y, X, multi=True)
    mgwr_bw = mgwr_selector.search(multi_bw_min=[2])
    mgwr_results = MGWR(coords, y, X, mgwr_selector).fit()

    return data, np.asarray(gwr_bw), gwr_results, np.asarray(mgwr_bw), mgwr_results


def diagnostics(result) -> dict[str, float]:
    return {
        "rss": _scalar(result.RSS),
        "enp": _scalar(result.ENP),
        "aic": float(result.aic),
        "aicc": float(result.aicc),
        "bic": float(result.bic),
        "r2": float(result.R2),
        "adj_r2": float(result.adj_R2),
    }


def build_comparison(actual: dict) -> dict:
    comparison = {"gwr": {}, "mgwr": {}}
    comparison["gwr"]["bandwidth_delta"] = (
        actual["gwr"]["bandwidth"] - REFERENCE["gwr"]["bandwidth"]
    )
    comparison["mgwr"]["bandwidth_deltas"] = [
        float(a - b)
        for a, b in zip(actual["mgwr"]["bandwidths"], REFERENCE["mgwr"]["bandwidths"])
    ]
    for model in ("gwr", "mgwr"):
        for key in ("rss", "enp", "aic", "aicc", "bic", "r2", "adj_r2"):
            comparison[model][f"{key}_delta"] = (
                actual[model][key] - REFERENCE[model][key]
            )
    return comparison


def save_outputs(data, gwr_bw, gwr_results, mgwr_bw, mgwr_results):
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    gwr_params = pd.DataFrame(gwr_results.params, columns=VARIABLE_NAMES)
    mgwr_params = pd.DataFrame(mgwr_results.params, columns=VARIABLE_NAMES)
    id_cols = data[["AreaKey", "ID", "X", "Y"]].reset_index(drop=True)
    pd.concat([id_cols, gwr_params], axis=1).to_csv(
        RESULT_DIR / "gwr_parameters.csv", index=False
    )
    pd.concat([id_cols, mgwr_params], axis=1).to_csv(
        RESULT_DIR / "mgwr_parameters.csv", index=False
    )

    actual = {
        "n": int(len(data)),
        "variables": {
            "y": "PctBach",
            "x": ["PctFB", "PctBlack", "PctRural"],
            "coords": ["X", "Y"],
        },
        "preprocessing": "z-score y and X columns with ddof=0",
        "kernel": "adaptive bisquare",
        "criterion": "AICc",
        "versions": {
            package: metadata.version(package)
            for package in ["mgwr", "libpysal", "numpy", "scipy", "pandas"]
        },
        "gwr": {
            "bandwidth": _scalar(gwr_bw),
            **diagnostics(gwr_results),
            "parameter_summary": parameter_summary(gwr_results.params),
        },
        "mgwr": {
            "bandwidths": [float(v) for v in np.asarray(mgwr_bw).reshape(-1)],
            **diagnostics(mgwr_results),
            "parameter_summary": parameter_summary(mgwr_results.params),
        },
    }

    payload = {
        "reference": REFERENCE,
        "actual": actual,
        "difference_actual_minus_reference": build_comparison(actual),
    }
    (RESULT_DIR / "summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    lines = [
        "# Georgia MGWR 2019 reproduction",
        "",
        f"- observations: {actual['n']}",
        f"- mgwr version: {actual['versions']['mgwr']}",
        f"- GWR bandwidth: {actual['gwr']['bandwidth']}",
        f"- MGWR bandwidths: {actual['mgwr']['bandwidths']}",
        f"- GWR AICc: {actual['gwr']['aicc']:.6f}",
        f"- MGWR AICc: {actual['mgwr']['aicc']:.6f}",
        f"- GWR R2: {actual['gwr']['r2']:.6f}",
        f"- MGWR R2: {actual['mgwr']['r2']:.6f}",
        "",
        "The machine-readable comparison with the historical notebook is in `summary.json`.",
    ]
    (RESULT_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    return actual


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vendor-data",
        action="store_true",
        help="copy canonical libpysal Georgia CSV/shapefile files into data/raw/georgia",
    )
    args = parser.parse_args()

    data_path = locate_data(args.vendor_data)
    data, gwr_bw, gwr_results, mgwr_bw, mgwr_results = fit_models(data_path)
    actual = save_outputs(data, gwr_bw, gwr_results, mgwr_bw, mgwr_results)

    print("Georgia reproduction complete")
    print("GWR bandwidth:", actual["gwr"]["bandwidth"])
    print("MGWR bandwidths:", actual["mgwr"]["bandwidths"])
    print("GWR AICc:", actual["gwr"]["aicc"])
    print("MGWR AICc:", actual["mgwr"]["aicc"])
    print("GWR R2:", actual["gwr"]["r2"])
    print("MGWR R2:", actual["mgwr"]["r2"])


if __name__ == "__main__":
    main()
