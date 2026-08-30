"""Minimal manual smoke run for the research baseline."""

import numpy as np

from georegime_gwr import BasicGWR, GRGWRBaseline


rng = np.random.default_rng(42)
n = 120
coords = rng.uniform(0, 10, size=(n, 2))
X = rng.normal(size=(n, 2))
regime = (coords[:, 0] >= 5.0).astype(int)
beta1 = np.where(regime == 0, 2.0, -2.0)
y = 1.0 + beta1 * X[:, 0] + 0.8 * X[:, 1] + rng.normal(0, 0.08, n)

base = BasicGWR(bandwidth=25).fit(X, y, coords)
model = GRGWRBaseline(n_regimes=2, bandwidth=25, random_state=42).fit(X, y, coords)

print("GWR RMSE:", np.sqrt(np.mean(base.residuals_**2)))
print("GR-GWR RMSE:", np.sqrt(np.mean(model.residuals_**2)))
print("GR-GWR regimes:", np.bincount(model.regimes_))
print("Objective history:", model.objective_history_)
