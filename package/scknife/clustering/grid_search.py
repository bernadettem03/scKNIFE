import re

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA as sklearn_PCA


def extract_h_graph(adata):
    """Extract numerically ordered ``H_graph_*`` factors into ``X_H_raw``.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData containing per-cell factors in ``adata.obs`` columns named
        ``H_graph_0``, ``H_graph_1``, and so on.

    Returns
    -------
    list of str
        The factor columns in numeric component order.

    Notes
    -----
    This function mutates ``adata`` by storing a cell-by-factor copy in
    ``adata.obsm['X_H_raw']``.
    """
    pattern = re.compile(r"^H_graph_(\d+)$")
    indexed = [(int(match.group(1)), column) for column in adata.obs.columns
               if (match := pattern.match(str(column)))]
    indexed.sort(key=lambda item: item[0])
    h_cols = [column for _, column in indexed]
    if not h_cols:
        raise ValueError("No H_graph columns found in adata.obs")
    component_ids = [index for index, _ in indexed]
    if len(component_ids) != len(set(component_ids)):
        raise ValueError("Duplicate numeric H_graph component identifiers found")
    H = adata.obs[h_cols].to_numpy(dtype=float)
    if not np.isfinite(H).all():
        raise ValueError("H_graph factors contain non-finite values")
    adata.obsm["X_H_raw"] = H.copy()
    return h_cols


def run_grid_search(adata, pca_dims, neighbors_list, resolutions, metric="cosine",
                    stage_label="Stage1", use_h_graph=False, random_state=42,
                    sample_size=2000):
    """Evaluate PCA, neighbor-count, and Leiden-resolution combinations.

    When ``use_h_graph=True``, PCA is computed from ``adata.obsm['X_H_raw']``.
    Each configuration builds a neighbor graph, runs Leiden, and evaluates
    cosine (or ``metric``) silhouette in the corresponding PCA space.

    Parameters
    ----------
    adata : anndata.AnnData
        Input object. ``X_H_raw`` must already exist when ``use_h_graph=True``.
    pca_dims : sequence of int
        PCA dimensions to evaluate.
    neighbors_list : sequence of int
        Neighbor counts to evaluate.
    resolutions : sequence of float
        Leiden resolutions to evaluate.
    metric : str, default="cosine"
        Distance metric for neighbor construction and silhouette evaluation.
    stage_label : str, default="Stage1"
        Label recorded in the returned table.
    use_h_graph : bool, default=False
        Use ``X_H_raw`` rather than expression values as the PCA input.
    random_state : int, default=42
        Seed for PCA, Leiden, and silhouette subsampling.
    sample_size : int or None, default=2000
        Maximum cells used for each silhouette calculation; ``None`` uses all.

    Returns
    -------
    pandas.DataFrame
        One row per configuration with parameters, cluster count, silhouette,
        generated label key, and stage label.

    Notes
    -----
    The function mutates ``adata`` by adding PCA representations, neighbor
    graphs, and one Leiden label column per evaluated configuration.
    """
    results = []
    total = len(pca_dims) * len(neighbors_list) * len(resolutions)
    cfg = 0

    print("=" * 100)
    print(f"GRID SEARCH {stage_label}: {total} configs")
    print("=" * 100)
    print(f"{'Config':<8} {'PCA':<6} {'N_neigh':<8} {'Resol':<8} {'Clusters':<10} {'Silhouette':<12}")
    print("-" * 100)

    for p in pca_dims:
        if use_h_graph:
            if "X_H_raw" not in adata.obsm:
                raise ValueError("X_H_raw missing. Call extract_h_graph() first.")
            X_pca = sklearn_PCA(n_components=p, svd_solver="auto",
                                random_state=random_state).fit_transform(adata.obsm["X_H_raw"])
            adata.obsm["X_pca"] = X_pca
        else:
            sc.pp.pca(adata, n_comps=p, use_highly_variable=False, svd_solver="arpack")
            X_pca = adata.obsm["X_pca"]

        for n in neighbors_list:
            neighbors_key = f"neighbors_pca{p}_n{n}"
            sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=n, metric=metric, key_added=neighbors_key)

            for r in resolutions:
                cfg += 1
                cluster_key = f"leiden_pca{p}_n{n}_r{r}"
                sc.tl.leiden(adata, resolution=r, key_added=cluster_key,
                             neighbors_key=neighbors_key, random_state=random_state)
                labels = adata.obs[cluster_key].astype(int).values
                ncl = len(np.unique(labels))

                n_sample = min(sample_size, adata.n_obs) if sample_size else None
                sil = silhouette_score(X_pca, labels, metric=metric,
                                       sample_size=n_sample, random_state=random_state)

                results.append({
                    "pca_dims": p,
                    "n_neighbors": int(n),
                    "resolution": float(r),
                    "n_clusters": int(ncl),
                    "silhouette": float(sil),
                    "cluster_key": cluster_key,
                    "stage": stage_label,
                })
                print(f"{cfg:<8} {p:<6} {int(n):<8} {r:<8.2f} {ncl:<10} {sil:<12.4f}")

    return pd.DataFrame(results)


def get_adjacent_n_neighbors(n, all_vals=(5, 10, 15)):
    idx = all_vals.index(n)
    adj = []
    if idx > 0:
        adj.append(all_vals[idx - 1])
    if idx < len(all_vals) - 1:
        adj.append(all_vals[idx + 1])
    return adj
