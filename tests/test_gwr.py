import numpy as np

from georegime_gwr import BasicGWR


def test_basic_gwr_runs_and_returns_one_parameter_vector_per_location():
    rng = np.random.default_rng(0)
    n = 40
    coords = rng.uniform(0, 10, size=(n, 2))
    X = rng.normal(size=(n, 2))
    y = 2.0 + 1.5 * X[:, 0] - 0.7 * X[:, 1] + rng.normal(0, 0.05, n)

    model = BasicGWR(bandwidth=15, kernel="bisquare").fit(X, y, coords)

    assert model.parameters_.shape == (n, 3)
    assert model.fitted_values_.shape == (n,)
    assert model.hat_matrix_.shape == (n, n)
    assert np.all(np.isfinite(model.parameters_))
    assert model.adaptive_ is True
    assert model.boundary_policy_ == "pygwrx"


def test_auto_policy_defaults_are_explicit():
    adaptive_default = BasicGWR(bandwidth="auto")
    assert adaptive_default._resolve_adaptive() is True
    assert adaptive_default._resolve_search_strategy(True) == "exhaustive"

    adaptive_compat = BasicGWR(
        bandwidth="auto",
        adaptive=True,
        search_strategy="mgwr_golden",
    )
    assert adaptive_compat._resolve_search_strategy(True) == "mgwr_golden"

    fixed_default = BasicGWR(bandwidth="auto", adaptive=False)
    assert fixed_default._resolve_adaptive() is False
    assert fixed_default._resolve_search_strategy(False) == "golden_section"
