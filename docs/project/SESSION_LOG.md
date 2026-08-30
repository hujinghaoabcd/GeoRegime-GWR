# SESSION_LOG

## 2026-08-30 — Repository initialization

- Created independent research repository `GeoRegime-GWR`.
- Added minimal `BasicGWR` implementation.
- Added readable `GRGWRBaseline` extracted from the current pyGWRx logic.
- Explicitly marked the baseline as non-final and allowed future algorithm replacement.
- Added cross-conversation handoff (`00_PROJECT_HANDOFF.md`).
- Added architecture index and current-status file.
- Added mathematical baseline specification and open research questions.
- Added staged experiment plan.
- Added initial ADR defining repository scope.
- Added lightweight tests, smoke script, and CI workflow.
- Recorded the current highest-priority design issue: whether the model should be formulated around an explicit spatial adjacency matrix `W` rather than hard-coded kNN+MST construction.

Future sessions should append a short dated entry whenever the model definition or experiment plan materially changes.
