"""Validate and extract the genes-by-cells fitting matrix."""
import numpy as np
from scipy import sparse


def extract_expression_matrix(adata, G, max_dense_bytes=2_000_000_000):
    """Retain exact graph gene names and materialize a bounded dense matrix.

    Input must already be normalized on a nonnegative scale. The limit covers
    the input array only; fitting requires additional working memory.
    """
    if not adata.var_names.is_unique or not adata.obs_names.is_unique:
        raise ValueError("Gene and cell names must be unique")
    valid = adata.var_names.isin(G.nodes)
    adata = adata[:, valid].copy()
    if not adata.n_obs or not adata.n_vars:
        raise ValueError("Need cells and at least one exact gene-name match in the graph")
    required = adata.n_obs * adata.n_vars * 8
    if required > max_dense_bytes:
        raise ValueError(f"Dense expression requires {required:,} bytes; subset genes/cells or explicitly increase max_dense_bytes")
    X = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    X = np.asarray(X, dtype=float).T.copy()
    if not np.isfinite(X).all() or (X < 0).any() or not (X > 0).any():
        raise ValueError("Expression must be finite, nonnegative, and not entirely zero")
    return X, adata
