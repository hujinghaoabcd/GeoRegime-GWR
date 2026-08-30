"""Reproduce the canonical Georgia standard-GWR example with the ``mgwr`` package.

This experiment intentionally uses only ``mgwr.gwr.GWR``.  MGWR is out of
scope at the current research stage.  The purpose is to establish a trusted
standard-GWR baseline before validating or modifying GR-GWR.

Specification
-------------
Dependent variable: PctBach
Explanatory variables: PctFB, PctBlack, PctRural
Coordinates: projected X, Y
Pre-processing: z-score y and each X column with ddof=0
Kernel: adaptive bisquare
Bandwidth criterion: AICc
Historical reference bandwidth: 117
"""

from __future__ import annotations

import argparse
import json
import shutil
from importlib import metadata
from pathlib import Path

import libpysal as ps
import numpy as np
import pandas as pd
from mgwr.gwr import GWR
from mgwr.sel_bw import Sel_BW

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = ROOT / "data" / "raw" / "georgia"
RESULT_DIR = ROOT / "results" / "reproduction" / "georgia_gwr_baseline"
VARIABLE_NAMES = ["Intercept", "PctFB", "PctBlack", "PctRural"]

REFERENCE = {
    "source": "canonical PySAL/mgwr Georgia GWR example",
    "n": 159,
    "bandwidth": 117.0,
    "rss": 51.186,
    "enp": 11.805,
    "aic": 296.616,
    "aicc": 299.051,
    "bic": 335.913,
    "r2": 0.678,
    "adj_r2": 0.652,
}


def _scalar(value) -> float:
    return float(np.asarray(value).reshape(-1)[0])


def vendor_georgia_data() -> Path:
    """Copy the canonical libpysal Georgia data into this research repository."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    csv_src = Path(ps.examples.get_path("GData_utm.csv"))
    csv_dst = DATA_DIR / "GData_utm.csv"
    shutil.copy2(csv_src, csv_dst)

    # Keep county polygons for later Queen/Rook spatial-weight experiments.
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
    out = {}
    for j, name in enumerate(VARIABLE_NAMES):
        v = np.asarray(params[:, j], dtype=float)
        out[name] = {
            "mean": float(v.mean()),
            "std": float(v.std(ddof=0)),
            "min": float(v.min()),
            "median": float(np.median(v)),
            "max": float(v.max()),
        }
    return out


def fit_gwr(data_path: Path):
    data = pd.read_csv(data_path)
    y = data["PctBach"].to_numpy(dtype=float).reshape((-1, 1))
    X = data[["PctFB", "PctBlack", "PctRural"]].to_numpy(dtype=float)
    coords = list(zip(data["X"].to_numpy(dtype=float), data["Y"].to_numpy(dtype=float)))

    # Match the canonical example exactly.
    X = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    y = (y - y.mean(axis=0)) / y.std(axis=0, ddof=0)

    selector = Sel_BW(coords, y, X)
    bandwidth = selector.search(bw_min=2)
    result = GWR(
        coords,
        y,
        X,
        bandwidth,
        fixed=False,
        kernel="bisquare",
        constant=True,
    ).fit()
    return data, float(bandwidth), result


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


def save_outputs(data: pd.DataFrame, bandwidth: float, result) -> dict:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    params = pd.DataFrame(result.params, columns=VARIABLE_NAMES)
    id_cols = data[["AreaKey", "ID", "X", "Y"]].reset_index(drop=True)
    pd.concat([id_cols, params], axis=1).to_csv(
        RESULT_DIR / "gwr_parameters.csv", index=False
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
        "bandwidth": bandwidth,
        **diagnostics(result),
        "parameter_summary": parameter_summary(result.params),
    }

    comparison = {
        key: actual[key] - REFERENCE[key]
        for key in ["bandwidth", "rss", "enp", "aic", "aicc", "bic", "r2", "adj_r2"]
    }
    payload = {
        "reference": REFERENCE,
        "actual": actual,
        "difference_actual_minus_reference": comparison,
    }
    (RESULT_DIR / "summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    lines = [
        "# Georgia standard GWR baseline",
        "",
        "This result uses `mgwr.gwr.GWR` only. MGWR is intentionally not run.",
        "",
        f"- observations: {actual['n']}",
        f"- mgwr package version: {actual['versions']['mgwr']}",
        f"- bandwidth: {actual['bandwidth']}",
        f"- RSS: {actual['rss']:.6f}",
        f"- ENP: {actual['enp']:.6f}",
        f"- AICc: {actual['aicc']:.6f}",
        f"- R2: {actual['r2']:.6f}",
        f"- adjusted R2: {actual['adj_r2']:.6f}",
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

    data, bandwidth, result = fit_gwr(locate_data(args.vendor_data))
    actual = save_outputs(data, bandwidth, result)
    print("Georgia standard-GWR reproduction complete")
    print("Bandwidth:", actual["bandwidth"])
    print("AICc:", actual["aicc"])
    print("R2:", actual["r2"])


if __name__ == "__main__":
    main()
