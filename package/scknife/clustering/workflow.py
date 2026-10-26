"""User-facing H-factor clustering workflow."""
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from .grid_search import extract_h_graph, run_grid_search


def _validate_fixed_parameters(adata, pca_dims, n_neighbors, resolution):
    if "X_H_raw" not in adata.obsm:
        raise ValueError("X_H_raw is missing; call extract_h_graph() first")
    n_cells, n_factors = adata.obsm["X_H_raw"].shape
    if not isinstance(pca_dims, (int, np.integer)) or not 1 <= pca_dims <= min(n_cells, n_factors):
        raise ValueError(f"pca_dims must be between 1 and {min(n_cells, n_factors)}")
    if not isinstance(n_neighbors, (int, np.integer)) or not 2 <= n_neighbors < n_cells:
        raise ValueError(f"n_neighbors must be between 2 and {n_cells - 1}")
    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError("resolution must be a positive finite number")


def _fit_h_graph_clustering(adata, *, pca_dims, n_neighbors, resolution,
                            metric, random_state, sample_size=2000):
    """Fit one validated H-PCA-neighbors-Leiden configuration."""
    _validate_fixed_parameters(adata, pca_dims, n_neighbors, resolution)
    adata.obsm["X_pca"] = PCA(
        n_components=int(pca_dims), svd_solver="auto", random_state=random_state,
    ).fit_transform(adata.obsm["X_H_raw"])
    sc.pp.neighbors(
        adata, use_rep="X_pca", n_neighbors=int(n_neighbors), metric=metric,
        key_added="nb_final",
    )
    sc.tl.leiden(
        adata, resolution=float(resolution), key_added="leiden_best",
        neighbors_key="nb_final", random_state=random_state,
    )
    labels = adata.obs["leiden_best"].astype(int).to_numpy()
    adata.obs["leiden_clusters"] = labels
    adata.obs["leiden_clusters_str"] = adata.obs["leiden_clusters"].astype(str)
    sc.tl.umap(adata, neighbors_key="nb_final", random_state=random_state)
    n_sample = min(sample_size, adata.n_obs) if sample_size else None
    silhouette = -1.0
    n_clusters = int(np.unique(labels).size)
    if 1 < n_clusters < adata.n_obs:
        silhouette = float(silhouette_score(
            adata.obsm["X_pca"], labels, metric=metric,
            sample_size=n_sample, random_state=random_state,
        ))
    return pd.DataFrame([{
        "n_clusters": n_clusters,
        "best_pca_dims": int(pca_dims),
        "best_n_neighbors": int(n_neighbors),
        "best_resolution": float(resolution),
        "metric": metric,
        "random_state": int(random_state),
        "silhouette": silhouette,
    }])


def _write_fixed_outputs(adata, summary, output_dir):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out / "clustering_summary.csv", index=False)
    adata.obs[["leiden_clusters"]].to_csv(out / "cell_clustering_assignment.csv")
    adata.write_h5ad(out / "clustered_with_H.h5ad")


def run_h_graph_clustering_fixed(adata, *, pca_dims, n_neighbors, resolution,
                                 metric="cosine", random_state=42,
                                 output_dir=None):
    """Run one fixed H-factor PCA, neighbor, Leiden, and UMAP configuration.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData containing ``H_graph_<integer>`` columns in ``adata.obs``.
    pca_dims : int
        Number of principal components computed from the H-factor matrix.
    n_neighbors : int
        Number of neighbors used to construct the graph in PCA space.
    resolution : float
        Positive Leiden resolution.
    metric : str, default="cosine"
        Distance metric used for neighbors and silhouette evaluation.
    random_state : int, default=42
        Seed used by PCA, Leiden, UMAP, and silhouette subsampling.
    output_dir : path-like, optional
        If supplied, save the summary, cell assignments, and clustered h5ad.

    Returns
    -------
    adata : anndata.AnnData
        The input object, mutated with ``X_H_raw``, ``X_pca``, ``X_umap``,
        ``leiden_best``, ``leiden_clusters``, and ``leiden_clusters_str``.
    summary : pandas.DataFrame
        One-row table containing the parameters, metric, seed, number of
        clusters, and silhouette score.
    """
    extract_h_graph(adata)
    summary = _fit_h_graph_clustering(
        adata, pca_dims=pca_dims, n_neighbors=n_neighbors,
        resolution=resolution, metric=metric, random_state=random_state,
    )
    if output_dir is not None:
        _write_fixed_outputs(adata, summary, output_dir)
    return adata, summary


def run_h_graph_clustering(adata, pca_dims=(10,15), neighbors=(5,10,15),
                           resolutions=(0.1,0.3,0.5,1.0), metric="cosine",
                           random_state=42, output_dir=None):
    """Select an H-space Leiden solution by maximum silhouette and compute UMAP.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData containing ``H_graph_<integer>`` columns in ``adata.obs``.
    pca_dims : sequence of int, default=(10, 15)
        PCA dimensions included in the grid.
    neighbors : sequence of int, default=(5, 10, 15)
        Neighbor counts included in the grid.
    resolutions : sequence of float, default=(0.1, 0.3, 0.5, 1.0)
        Leiden resolutions included in the grid.
    metric : str, default="cosine"
        Distance metric used for neighbors and silhouette evaluation.
    random_state : int, default=42
        Seed used by PCA, Leiden, UMAP, and silhouette subsampling.
    output_dir : path-like, optional
        If supplied, save the grid, selected summary, assignments, and h5ad.

    Returns
    -------
    adata : anndata.AnnData
        Mutated input containing the final PCA, neighbor graph, labels, and UMAP.
    evaluation : pandas.DataFrame
        Complete grid with one row per evaluated configuration.
    summary : pandas.DataFrame
        One-row table describing the maximum-silhouette final configuration.

    Notes
    -----
    This is a label-free maximum-silhouette selection workflow. Use
    :func:`run_h_graph_clustering_fixed` when parameters are already known.
    """
    extract_h_graph(adata)
    evaluation=run_grid_search(adata,list(pca_dims),list(neighbors),list(resolutions),
                               metric=metric,stage_label="Unsupervised",use_h_graph=True,
                               random_state=random_state)
    best=evaluation.loc[evaluation["silhouette"].idxmax()]
    summary = _fit_h_graph_clustering(
        adata, pca_dims=int(best.pca_dims), n_neighbors=int(best.n_neighbors),
        resolution=float(best.resolution), metric=metric,
        random_state=random_state,
    )
    if output_dir is not None:
        out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)
        summary.to_csv(out/'best_clustering.csv',index=False)
        evaluation.to_csv(out/'clustering_evaluation.csv',index=False)
        evaluation.sort_values('silhouette',ascending=False).head(10).to_csv(out/'clustering_top10_summary.csv',index=False)
        adata.obs[["leiden_clusters"]].to_csv(out/'cell_clustering_assignment.csv')
        adata.write_h5ad(out/'clustered_with_H.h5ad')
    return adata,evaluation,summary
