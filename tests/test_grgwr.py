import numpy as np

from georegime_gwr import GRGWRBaseline


def _two_regime_data(seed=1, n=120):
    rng = np.random.default_rng(seed)
    coords = rng.uniform(0, 10, size=(n, 2))
    X = rng.normal(size=(n, 2))
    true_regime = (coords[:, 0] >= 5.0).astype(int)
    beta1 = np.where(true_regime == 0, 2.0, -2.0)
    y = 1.0 + beta1 * X[:, 0] + 0.8 * X[:, 1] + rng.normal(0, 0.08, n)
    return X, y, coords, true_regime


def test_grgwr_baseline_runs_without_refinement():
    X, y, coords, _ = _two_regime_data()
    model = GRGWRBaseline(
        n_regimes=2,
        bandwidth=25,
        n_neighbors=8,
        refine=False,
        random_state=42,
    ).fit(X, y, coords)

    assert model.regimes_.shape == (X.shape[0],)
    assert model.parameters_.shape == (X.shape[0], X.shape[1] + 1)
    assert len(np.unique(model.regimes_)) == 2
    assert np.all(np.isfinite(model.fitted_values_))


def test_accepted_refinement_objective_is_nonincreasing():
    X, y, coords, _ = _two_regime_data(seed=2)
    model = GRGWRBaseline(
        n_regimes=2,
        bandwidth=25,
        n_neighbors=8,
        refine=True,
        max_iter=5,
        random_state=42,
    ).fit(X, y, coords)

    history = np.asarray(model.objective_history_)
    assert np.all(np.diff(history) <= 1e-8)


def test_each_location_still_has_its_own_local_parameter_vector():
    X, y, coords, _ = _two_regime_data(seed=3)
    model = GRGWRBaseline(
        n_regimes=2,
        bandwidth=25,
        refine=False,
    ).fit(X, y, coords)

    assert model.parameters_.shape[0] == X.shape[0]
    # GR-GWR is not "one coefficient vector per regime".
    assert np.unique(np.round(model.parameters_[:, 1], 6)).size > 2
