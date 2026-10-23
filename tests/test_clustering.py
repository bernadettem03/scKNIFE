import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scknife.clustering import (
    extract_h_graph,
    run_h_graph_clustering,
    run_h_graph_clustering_fixed,
)


def make_h_adata(n_cells=48):
    rng = np.random.RandomState(7)
    half = n_cells // 2
    first = rng.normal(loc=-2.0, scale=0.25, size=(half, 4))
    second = rng.normal(loc=2.0, scale=0.25, size=(n_cells - half, 4))
    factors = np.vstack([first, second])
    obs = pd.DataFrame(
        {
            "H_graph_2": factors[:, 2],
            "H_graph_0": factors[:, 0],
            "H_graph_3": factors[:, 3],
            "H_graph_1": factors[:, 1],
        },
        index=[f"cell_{index}" for index in range(n_cells)],
    )
    return AnnData(X=np.zeros((n_cells, 1), dtype=np.float32), obs=obs), factors


def test_extract_h_graph_uses_numeric_component_order():
    adata, factors = make_h_adata()
    columns = extract_h_graph(adata)
    assert columns == ["H_graph_0", "H_graph_1", "H_graph_2", "H_graph_3"]
    np.testing.assert_allclose(adata.obsm["X_H_raw"], factors)


def test_fixed_clustering_populates_outputs_and_is_reproducible(tmp_path):
    first, _ = make_h_adata()
    second, _ = make_h_adata()
    first, summary = run_h_graph_clustering_fixed(
        first, pca_dims=2, n_neighbors=5, resolution=0.3,
        random_state=42, output_dir=tmp_path,
    )
    second, second_summary = run_h_graph_clustering_fixed(
        second, pca_dims=2, n_neighbors=5, resolution=0.3,
        random_state=42,
    )

    assert {"X_H_raw", "X_pca", "X_umap"}.issubset(first.obsm)
    assert {"leiden_best", "leiden_clusters", "leiden_clusters_str"}.issubset(first.obs)
    assert "nb_final" in first.uns
    assert summary.loc[0, "best_pca_dims"] == 2
    assert summary.loc[0, "best_n_neighbors"] == 5
    assert summary.loc[0, "best_resolution"] == pytest.approx(0.3)
    assert summary.loc[0, "metric"] == "cosine"
    assert summary.loc[0, "random_state"] == 42
    np.testing.assert_array_equal(
        first.obs["leiden_clusters"], second.obs["leiden_clusters"]
    )
    pd.testing.assert_frame_equal(summary, second_summary)
    for filename in (
        "clustering_summary.csv",
        "cell_clustering_assignment.csv",
        "clustered_with_H.h5ad",
    ):
        assert (tmp_path / filename).is_file()


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"pca_dims": 0, "n_neighbors": 5, "resolution": 0.3}, "pca_dims"),
        ({"pca_dims": 2, "n_neighbors": 1, "resolution": 0.3}, "n_neighbors"),
        ({"pca_dims": 2, "n_neighbors": 5, "resolution": 0}, "resolution"),
    ],
)
def test_fixed_clustering_validates_parameters(kwargs, message):
    adata, _ = make_h_adata()
    with pytest.raises(ValueError, match=message):
        run_h_graph_clustering_fixed(adata, **kwargs)


def test_grid_workflow_refits_selected_configuration(monkeypatch):
    adata, _ = make_h_adata()
    selected = pd.DataFrame(
        [
            {"pca_dims": 2, "n_neighbors": 5, "resolution": 0.2,
             "n_clusters": 2, "silhouette": 0.4},
            {"pca_dims": 3, "n_neighbors": 7, "resolution": 0.5,
             "n_clusters": 2, "silhouette": 0.8},
        ]
    )
    captured = {}

    def fake_grid(*args, **kwargs):
        return selected

    def fake_fit(adata_arg, **kwargs):
        captured.update(kwargs)
        return pd.DataFrame([{"n_clusters": 2, "best_pca_dims": kwargs["pca_dims"],
                              "best_n_neighbors": kwargs["n_neighbors"],
                              "best_resolution": kwargs["resolution"],
                              "metric": kwargs["metric"],
                              "random_state": kwargs["random_state"],
                              "silhouette": 0.8}])

    monkeypatch.setattr("scknife.clustering.workflow.run_grid_search", fake_grid)
    monkeypatch.setattr("scknife.clustering.workflow._fit_h_graph_clustering", fake_fit)
    _, evaluation, summary = run_h_graph_clustering(
        adata, pca_dims=(2, 3), neighbors=(5, 7), resolutions=(0.2, 0.5),
        random_state=11,
    )

    assert evaluation is selected
    assert captured == {
        "pca_dims": 3,
        "n_neighbors": 7,
        "resolution": 0.5,
        "metric": "cosine",
        "random_state": 11,
    }
    assert summary.loc[0, "best_pca_dims"] == 3
