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
