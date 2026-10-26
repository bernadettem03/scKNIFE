"""
Advanced tuning criteria for scKNIFE hyperparameter optimization.

Provides:
  - reconstruction_correlation : Pearson corr(X, X_hat)
  - leiden_silhouette          : best silhouette over Leiden resolutions
  - evaluate_h_graph_clustering: direct H_graph clustering (no PCA) with silhouette
  - factor_stability           : factor matching + gene/go/hallmark Jaccard
  - find_elbow_kneedle         : Kneedle algorithm for reconstruction elbow
  - compute_composite_score    : equal-weight normalized composite
  - write_tuning_report        : full text report
"""

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score

from scknife.evaluation.reconstruction import evaluate_reconstruction


STABILITY_CORR_THRESH = 0.7

LEIDEN_RESOLUTIONS = [0.5, 1.0, 1.2]


def reconstruction_correlation(X, Rg, W, H):
    """Criterion 1: Pearson correlation + errors."""
    corr, fro, rel = evaluate_reconstruction(X, Rg, W, H)
    return corr, fro, rel


def leiden_silhouette(H_mat, resolutions=None, random_state=0):
    import scanpy as sc
    import anndata

    if resolutions is None:
        resolutions = LEIDEN_RESOLUTIONS

    n_cells, k = H_mat.shape
    if n_cells < 5 or k < 2:
        return -1.0, 1.0, 1, None

    adata_tmp = anndata.AnnData(H_mat)
    try:
        sc.pp.neighbors(adata_tmp, n_neighbors=min(15, n_cells - 1), use_rep='X')
    except Exception:
        return -1.0, 1.0, 1, None

    best_sil = -1.0
    best_res = 1.0
    best_n = 1
    best_labels = None

    for res in resolutions:
        try:
            sc.tl.leiden(adata_tmp, resolution=res, random_state=random_state)
            labels = adata_tmp.obs['leiden'].astype(int).values
            n_unique = len(np.unique(labels))
            if n_unique > 1 and n_unique < n_cells:
                sil = silhouette_score(H_mat, labels)
                if sil > best_sil:
                    best_sil = sil
                    best_res = res
                    best_n = n_unique
                    best_labels = labels
        except Exception:
            continue

    return best_sil, best_res, best_n, best_labels


def evaluate_h_graph_clustering(H_mat, n_neighbors=20, metric="cosine",
                                resolution=0.5, random_state=42, sample_size=2000):
    """
    Cluster cells using H_graph features directly (no PCA) and evaluate via silhouette.

    Steps:
      1. Build neighbor graph on H_mat with specified metric.
      2. Run Leiden clustering at the given resolution.
      3. Compute silhouette score on H_mat with same metric.

    Parameters
    ----------
    H_mat : np.ndarray
        H_graph matrix of shape (n_cells, n_components).
    n_neighbors : int
        Number of neighbors for graph construction (default 20).
    metric : str
        Distance metric for neighbor graph and silhouette (default "cosine").
    resolution : float
        Leiden resolution (default 0.5).
    random_state : int
        Random seed for Leiden (default 42).
    sample_size : int or None
        Cells to sample for silhouette (default 2000). None = all.

    Returns
    -------
    dict with keys:
        "silhouette"  : float, best silhouette (-1.0 if <2 clusters)
        "n_clusters"  : int, number of clusters found
        "resolution"  : float, the resolution used
    """
    import scanpy as sc
    import anndata

    n_cells, k = H_mat.shape
    if n_cells < 5 or k < 2:
        return {"silhouette": -1.0, "n_clusters": 1, "resolution": resolution}

    adata = anndata.AnnData(X=H_mat)
    adata.obsm["X_H_raw"] = H_mat
    try:
        sc.pp.neighbors(adata, use_rep="X_H_raw", n_neighbors=min(n_neighbors, n_cells - 1),
                        metric=metric)
    except Exception:
        return {"silhouette": -1.0, "n_clusters": 1, "resolution": resolution}

    try:
        sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
        labels = adata.obs["leiden"].astype(int).values
        n_unique = len(np.unique(labels))
        if n_unique < 2 or n_unique >= n_cells:
            return {"silhouette": -1.0, "n_clusters": n_unique, "resolution": resolution}
        ss = min(sample_size, n_cells) if sample_size else n_cells
        sil = silhouette_score(H_mat, labels, metric=metric, sample_size=ss)
    except Exception:
        return {"silhouette": -1.0, "n_clusters": 1, "resolution": resolution}

    return {"silhouette": float(sil), "n_clusters": int(n_unique), "resolution": resolution}


def _greedy_match(Wa, Wb, k, n_sample):
    """Greedy factor matching between two W matrices. Returns matched corrs + Jaccards."""
    if k > n_sample:
        idx = np.random.RandomState(0).choice(k, size=n_sample, replace=False)
        Wa_s, Wb_s = Wa[:, idx], Wb[:, idx]
        sample_k = n_sample
    else:
        Wa_s, Wb_s = Wa, Wb
        sample_k = k

    corr_matrix = np.corrcoef(Wa_s.T, Wb_s.T)[:sample_k, sample_k:]
    used_b = set()
    matched_corrs = []
    for a_idx in range(sample_k):
        best_b = None
        best_c = -2
        for b_idx in range(sample_k):
            if b_idx in used_b:
                continue
            c = corr_matrix[a_idx, b_idx]
            if c > best_c:
                best_c = c
                best_b = b_idx
        if best_b is not None and not np.isnan(best_c):
            used_b.add(best_b)
            matched_corrs.append(best_c)

    return matched_corrs, Wa_s, Wb_s, list(used_b)


def _jaccard_top(W_a, W_b, matched_b_indices, top_n):
    """Jaccard of top-N genes for matched factors."""
    jaccs = []
    for a_idx in range(W_a.shape[1]):
        b_idx = matched_b_indices[a_idx]
        top_a = set(np.argsort(W_a[:, a_idx])[-top_n:])
        top_b = set(np.argsort(W_b[:, b_idx])[-top_n:])
        intersection = len(top_a & top_b)
        union = len(top_a | top_b)
        jaccs.append(intersection / union if union > 0 else 0.0)
    return float(np.mean(jaccs)) if jaccs else 0.0


def _subset_jaccard(W_a, W_b, matched_b_indices, mask, top_n):
    """Jaccard of top-N entries within a masked subset."""
    jaccs = []
    for a_idx in range(W_a.shape[1]):
        b_idx = matched_b_indices[a_idx]
        vals_a = np.where(mask, W_a[:, a_idx], 0)
        vals_b = np.where(mask, W_b[:, b_idx], 0)
        top_a = set(np.argsort(vals_a)[-top_n:])
        top_b = set(np.argsort(vals_b)[-top_n:])
        intersection = len(top_a & top_b)
        union = len(top_a | top_b)
        jaccs.append(intersection / union if union > 0 else 0.0)
    return float(np.mean(jaccs)) if jaccs else 0.0


def factor_stability(W_list, node_kinds=None, n_sample=10, top_n=50):
    """Criterion 3: multi-level factor stability across seeds.

    For each pair of seeds:
      1. Match factors greedily by max Pearson correlation.
      2. Compute:
         - factor_corr    : mean correlation of matched factors
         - gene_jacc      : Jaccard of top-N genes
         - go_jacc        : Jaccard of top-N GO nodes
         - hallmark_jacc  : Jaccard of top-N Hallmark nodes
         - combined       : mean of factor_corr + all Jaccards

    Returns dict with all metrics + combined stability score.
    """
    n_seeds = len(W_list)
    if n_seeds < 2:
        return {
            "factor_corr": np.nan, "gene_jacc": np.nan,
            "go_jacc": np.nan, "hallmark_jacc": np.nan,
            "combined_stability": np.nan,
        }

    n_nodes, k = W_list[0].shape

    # Precompute masks
    masks = {}
    if node_kinds is not None:
        kinds = list(node_kinds)
        masks['gene'] = np.array([v == 'gene' for v in kinds])
        masks['go'] = np.array([v == 'go' for v in kinds])
        masks['hallmark'] = np.array([v == 'hallmark' for v in kinds])

    corr_all = []
    gene_jacc_all = []
    go_jacc_all = []
    hall_jacc_all = []

    for i in range(n_seeds):
        for j in range(i + 1, n_seeds):
            matched_corrs, Wa_s, Wb_s, b_idx_list = _greedy_match(
                W_list[i], W_list[j], k, n_sample)

            if not matched_corrs:
                continue

            corr_all.extend(matched_corrs)
            matched_b = {a: b for a, b in enumerate(b_idx_list)}

            # Overall gene Jaccard (on full W)
            gene_jacc_all.append(_jaccard_top(Wa_s, Wb_s, matched_b, top_n))

            # Subset Jaccards
            if 'go' in masks and masks['go'].any():
                go_jacc_all.append(_subset_jaccard(Wa_s, Wb_s, matched_b, masks['go'], top_n))
            if 'hallmark' in masks and masks['hallmark'].any():
                hall_jacc_all.append(_subset_jaccard(Wa_s, Wb_s, matched_b, masks['hallmark'], top_n))

    mean_corr = float(np.mean(corr_all)) if corr_all else np.nan
    mean_gene = float(np.mean(gene_jacc_all)) if gene_jacc_all else np.nan
    mean_go = float(np.mean(go_jacc_all)) if go_jacc_all else np.nan
    mean_hall = float(np.mean(hall_jacc_all)) if hall_jacc_all else np.nan

    components = [v for v in [mean_corr, mean_gene, mean_go, mean_hall] if not np.isnan(v)]
    combined = float(np.mean(components)) if components else np.nan

    return {
        "factor_corr": mean_corr,
        "gene_jacc": mean_gene,
        "go_jacc": mean_go,
        "hallmark_jacc": mean_hall,
        "combined_stability": combined,
    }


def find_elbow_kneedle(x_vals, y_vals, curve="concave", sensitivity=1.0):
    """Simple Kneedle algorithm to find the elbow point.

    For a concave curve (reconstruction increases with k and plateaus),
    finds the point of maximum curvature.

    Returns
    -------
    elbow_x : float
        The x-value at the elbow.
    elbow_idx : int
        Index in the original arrays.
    """
    if len(x_vals) < 3:
        return x_vals[-1], len(x_vals) - 1

    x = np.array(x_vals, dtype=float)
    y = np.array(y_vals, dtype=float)

    # Normalize to [0, 1]
    x_norm = (x - x.min()) / (x.max() - x.min() + 1e-10)
    y_norm = (y - y.min()) / (y.max() - y.min() + 1e-10)

    # Difference between the curve and the line connecting endpoints
    line = np.linspace(y_norm[0], y_norm[-1], len(y_norm))
    diff = y_norm - line

    if curve == "convex":
        diff = -diff

    # Elbow is point with max distance (for concave: max positive diff)
    # For reconstruction correlation which plateaus, we want the point
    # where adding more k gives diminishing returns.
    # Use sensitivity to find earliest point with diff near max
    max_diff = np.max(diff)
    threshold = sensitivity * max_diff

    elbow_idx = 0
    for i in range(len(diff)):
        if diff[i] >= threshold:
            if i > 0:
                elbow_idx = i - 1
            else:
                elbow_idx = 0
            break
    else:
        elbow_idx = len(x) - 1

    return x[elbow_idx], elbow_idx


def compute_composite_score(summary_df, criteria_columns):
    """Normalize each criterion to [0,1] and compute equal-weight composite."""
    df = summary_df.copy()
    norm_cols = []
    for col in criteria_columns:
        vals = df[col].values.astype(float)
        vmin, vmax = np.nanmin(vals), np.nanmax(vals)
        if vmax - vmin > 1e-10:
            df[f"{col}_norm"] = (vals - vmin) / (vmax - vmin)
        else:
            df[f"{col}_norm"] = 0.5
        norm_cols.append(f"{col}_norm")

    df["composite"] = df[norm_cols].mean(axis=1)
    return df


def write_tuning_report(report_path, param_name, sweep_values, seed_list,
                        all_seed_results, summary_df, best_val,
                        selection_note="", elbow_info=None):
    """Write a human-readable text report."""
    lines = []
    lines.append("=" * 80)
    lines.append(f"scKNIFE TUNING REPORT — {param_name}")
    lines.append("=" * 80)
    if selection_note:
        lines.append(f"Selection: {selection_note}")
    if elbow_info:
        lines.append(f"Elbow: k={elbow_info.get('elbow_k')}  "
                      f"recon_corr={elbow_info.get('elbow_corr', 0):.4f}")
    lines.append(f"Values tested: {len(sweep_values)}")
    lines.append(f"Seeds per value: {len(seed_list)}")
    lines.append(f"Total runs: {len(sweep_values) * len(seed_list)}")
    lines.append("")

    for val in sweep_values:
        lines.append("-" * 60)
        lines.append(f"  {param_name} = {val}")
        lines.append("-" * 60)
        seed_rows = all_seed_results.get(val, [])
        for sr in seed_rows:
            lines.append(
                f"    seed={sr['seed']:2d}  "
                f"corr={sr.get('recon_corr', 0):.4f}  "
                f"rel_err={sr.get('rel_error', 0):.4f}  "
                f"sil={sr.get('silhouette', 0):.4f}  "
                f"f_corr={sr.get('stability_factor_corr', 0):.4f}  "
                f"g_jacc={sr.get('stability_gene_jacc', 0):.4f}  "
                f"go_jacc={sr.get('stability_go_jacc', 0):.4f}  "
                f"hall_jacc={sr.get('stability_hallmark_jacc', 0):.4f}"
            )
        lines.append("")

    lines.append("=" * 60)
    lines.append("  AGGREGATED SUMMARY")
    lines.append("=" * 60)
    for _, row in summary_df.iterrows():
        comp = row.get("composite", 0)
        vals = [f"{c}={row[c]:.4f}" for c in ["recon_corr", "silhouette", "combined_stability"]
                if c in row and not np.isnan(row.get(c, np.nan))]
        lines.append(f"    {str(row.get(param_name, '?')):>5s}  "
                      f"composite={comp:.4f}  "
                      f"{'  '.join(vals)}")
    lines.append("")

    lines.append(f"  >>> BEST {param_name} = {best_val} <<<")

    # Stability check
    if 'combined_stability' in summary_df.columns:
        stab = summary_df.loc[summary_df[param_name] == best_val, 'combined_stability'].values
        if len(stab) > 0 and not np.isnan(stab[0]):
            passed = stab[0] >= STABILITY_CORR_THRESH
            lines.append(f"  >>> Combined stability: {stab[0]:.4f} "
                          f"({'PASS' if passed else 'WARNING: BELOW THRESHOLD'}) <<<")
    lines.append("=" * 80)

    report_text = "\n".join(lines)
    with open(report_path, 'w') as f:
        f.write(report_text)
    return report_text
