import gc
import numpy as np
import pandas as pd
import anndata


def compute_activity_scores(adata, batch_size=10000):
    print("Computing activity scores from W and H...")
    W = adata.uns["W_graph"]
    import re
    if not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")
    H_cols = [col for col in adata.obs.columns if str(col).startswith("H_graph_")]
    if any(re.fullmatch(r"H_graph_\d+", str(col)) is None for col in H_cols):
        raise ValueError("Invalid H_graph factor column name")
    H_cols = sorted(H_cols, key=lambda col: int(col.rsplit("_", 1)[1]))
    if H_cols != [f"H_graph_{i}" for i in range(W.shape[1])]:
        raise ValueError("H_graph factors must occur exactly once, numbered 0 through W.shape[1]-1")
    H = adata.obs[H_cols].to_numpy(dtype=float).T
    if not np.isfinite(W).all() or not np.isfinite(H).all():
        raise ValueError("W and H must contain finite values")
    if len(adata.uns["node_names"]) != W.shape[0]:
        raise ValueError("node_names must match W rows")
    cell_names = adata.obs_names.tolist()
    node_names = adata.uns['node_names']

    n_cells = H.shape[1]
    n_nodes = W.shape[0]

    # Batch computation to avoid OOM on large datasets
    if n_cells > batch_size:
        n_batches = (n_cells + batch_size - 1) // batch_size
        activity_matrix = np.empty((n_nodes, n_cells), dtype=np.float32)
        for start in range(0, n_cells, batch_size):
            end = min(start + batch_size, n_cells)
            activity_matrix[:, start:end] = (W @ H[:, start:end]).astype(np.float32)
            print(f"  Activity batch {start//batch_size + 1}/{n_batches}: cells {start}-{end}")
            gc.collect()
    else:
        activity_matrix = W @ H

    activity_df = pd.DataFrame(activity_matrix, index=node_names, columns=cell_names)
    adata.uns['activity_scores'] = activity_df
    print("Stored activity scores in `adata.uns['activity_scores']`")
    return adata
