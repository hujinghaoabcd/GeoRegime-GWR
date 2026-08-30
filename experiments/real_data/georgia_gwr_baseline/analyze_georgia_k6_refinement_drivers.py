"""Diagnose why the final K=6 label-changed counties moved.

This is a diagnostic for the current exploratory refinement baseline.  It does
not change labels or refit the model.  For each county whose final label differs
from its initial label, we identify the last accepted move that put the county
into its final regime and decompose the local cost improvement into:

    LOO regression gain + lambda * Queen-boundary gain.

Positive ``loo_error_gain`` means the candidate regime predicted the county
better in leave-one-out local regression.  Positive ``boundary_gain`` means the
move reduced the number of Queen-neighbour label disagreements.
"""

from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "real_data" / "georgia_k6_unified_label_refinement"
MOVES = OUT / "accepted_label_moves.csv"
CHANGED = OUT / "changed_counties.csv"
LABEL_PATH = OUT / "label_path.csv"
LAMBDA_BOUNDARY = 1.0
TOL = 1e-12


def _classify(error_gain: float, boundary_gain: int) -> str:
    error_pos = error_gain > TOL
    boundary_pos = boundary_gain > 0
    boundary_neg = boundary_gain < 0
    if error_pos and boundary_pos:
        return "joint_regression_and_boundary"
    if error_pos and boundary_gain == 0:
        return "regression_only"
    if error_gain <= TOL and boundary_pos:
        return "boundary_overcame_worse_or_equal_LOO"
    if error_pos and boundary_neg:
        return "regression_overcame_boundary_cost"
    return "other"


def main() -> None:
    moves = pd.read_csv(MOVES)
    changed = pd.read_csv(CHANGED)
    path = pd.read_csv(LABEL_PATH)

    if len(changed) != 35:
        raise RuntimeError(f"Expected 35 final changed counties, got {len(changed)}")

    moves["AreaKey"] = moves["AreaKey"].astype(int)
    changed["AreaKey"] = changed["AreaKey"].astype(int)
    path["AreaKey"] = path["AreaKey"].astype(int)

    rows = []
    for c in changed.itertuples(index=False):
        area = int(c.AreaKey)
        county_moves = moves.loc[moves["AreaKey"] == area].sort_values("iteration")
        if county_moves.empty:
            raise RuntimeError(f"No accepted move found for final changed county {area}")

        decisive = county_moves.iloc[-1]
        if int(decisive["to_regime"]) != int(c.final_regime):
            raise RuntimeError(
                f"Last move for {area} ends in {int(decisive['to_regime'])}, "
                f"not final regime {int(c.final_regime)}"
            )

        error_gain = float(decisive["current_loo_error"] - decisive["best_loo_error"])
        boundary_gain = int(
            decisive["current_boundary_disagreement"]
            - decisive["best_boundary_disagreement"]
        )
        boundary_term_gain = LAMBDA_BOUNDARY * boundary_gain
        total_gain = float(decisive["current_cost"] - decisive["best_cost"])
        reconstructed = error_gain + boundary_term_gain
        if not np.isclose(total_gain, reconstructed, atol=1e-10, rtol=1e-10):
            raise RuntimeError(f"Cost decomposition mismatch for {area}")

        seq = [int(c.initial_regime)] + county_moves["to_regime"].astype(int).tolist()
        move_path = "→".join(map(str, seq))

        rows.append(
            {
                "AreaKey": area,
                "initial_regime": int(c.initial_regime),
                "final_regime": int(c.final_regime),
                "accepted_moves_for_county": int(len(county_moves)),
                "move_path": move_path,
                "decisive_iteration": int(decisive["iteration"]),
                "decisive_from_regime": int(decisive["from_regime"]),
                "decisive_to_regime": int(decisive["to_regime"]),
                "current_loo_error": float(decisive["current_loo_error"]),
                "candidate_loo_error": float(decisive["best_loo_error"]),
                "loo_error_gain": error_gain,
                "current_boundary_disagreement": int(decisive["current_boundary_disagreement"]),
                "candidate_boundary_disagreement": int(decisive["best_boundary_disagreement"]),
                "boundary_gain": boundary_gain,
                "boundary_term_gain_lambda1": boundary_term_gain,
                "total_local_cost_gain": total_gain,
                "regression_supported": bool(error_gain > TOL),
                "boundary_supported": bool(boundary_gain > 0),
                "driver_class": _classify(error_gain, boundary_gain),
            }
        )

    detail = pd.DataFrame(rows).sort_values(
        ["initial_regime", "final_regime", "AreaKey"]
    )
    detail.to_csv(OUT / "final_changed_driver_detail.csv", index=False)

    class_counts = detail["driver_class"].value_counts().to_dict()
    reg_supported = int(detail["regression_supported"].sum())
    bnd_supported = int(detail["boundary_supported"].sum())
    both_supported = int((detail["regression_supported"] & detail["boundary_supported"]).sum())

    r34 = detail[(detail["initial_regime"] == 3) & (detail["final_regime"] == 4)].copy()
    if len(r34) != 11:
        raise RuntimeError(f"Expected 11 final R3→R4 counties, got {len(r34)}")
    r34.to_csv(OUT / "r3_to_r4_driver_detail.csv", index=False)

    summary = {
        "status": "diagnostic_only_current_exploratory_refinement",
        "lambda_boundary": LAMBDA_BOUNDARY,
        "final_changed_counties": int(len(detail)),
        "decisive_move_driver_counts": {str(k): int(v) for k, v in class_counts.items()},
        "regression_supported_count": reg_supported,
        "boundary_supported_count": bnd_supported,
        "both_regression_and_boundary_supported_count": both_supported,
        "count_with_multiple_accepted_moves": int((detail["accepted_moves_for_county"] > 1).sum()),
        "mean_loo_error_gain": float(detail["loo_error_gain"].mean()),
        "median_loo_error_gain": float(detail["loo_error_gain"].median()),
        "mean_boundary_gain": float(detail["boundary_gain"].mean()),
        "mean_total_local_cost_gain": float(detail["total_local_cost_gain"].mean()),
        "r3_to_r4": {
            "count": int(len(r34)),
            "regression_supported_count": int(r34["regression_supported"].sum()),
            "boundary_supported_count": int(r34["boundary_supported"].sum()),
            "both_supported_count": int((r34["regression_supported"] & r34["boundary_supported"]).sum()),
            "driver_counts": {
                str(k): int(v) for k, v in r34["driver_class"].value_counts().to_dict().items()
            },
            "mean_loo_error_gain": float(r34["loo_error_gain"].mean()),
            "median_loo_error_gain": float(r34["loo_error_gain"].median()),
            "mean_boundary_gain": float(r34["boundary_gain"].mean()),
            "mean_total_local_cost_gain": float(r34["total_local_cost_gain"].mean()),
        },
        "interpretation_note": (
            "Driver labels refer to the decisive last accepted local move into the final regime. "
            "They do not prove causal geographic mechanisms and do not validate lambda=1.0."
        ),
    }
    (OUT / "refinement_driver_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))
    print("\nR3 -> R4 decisive-move detail:\n")
    print(r34.to_string(index=False))


if __name__ == "__main__":
    main()
