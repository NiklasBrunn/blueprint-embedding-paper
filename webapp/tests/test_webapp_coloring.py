import numpy as np
from webapp.coloring import scheme_for_obs, gene_log1p


def test_obs_scheme_categorical_vs_numeric():
    assert scheme_for_obs(["a", "b", "a"])["type"] == "categorical"
    assert scheme_for_obs([0.1, 0.2, 0.3])["type"] == "sequential"


def test_gene_log1p_matches_manual():
    X = np.array([[0.0, 1.0], [3.0, 0.0]])
    var_names = ["GeneA", "GeneB"]
    out = gene_log1p(X, var_names, "genea")  # case-insensitive
    assert np.allclose(out, np.log1p([0.0, 3.0]))
