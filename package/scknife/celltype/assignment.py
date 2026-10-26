"""Marker-activity helpers for assigning cell types to scKNIFE clusters."""

import pandas as pd


def compute_marker_anova_f(markers, activity_df, cluster_labels):
    """Rank marker activities by their separation across clusters."""
    k = cluster_labels.nunique()
    n = len(cluster_labels)
    grand_mean = activity_df.loc[markers].mean(axis=1)
    results = []
    for marker in markers:
        values = activity_df.loc[marker].values
        between_sum_squares = 0.0
        within_sum_squares = 0.0
        n_clusters_high = 0
        for cluster in cluster_labels.unique():
            mask = cluster_labels.values == cluster
            group = values[mask]
            group_mean = group.mean()
            between_sum_squares += len(group) * (group_mean - grand_mean[marker]) ** 2
            within_sum_squares += ((group - group_mean) ** 2).sum()
            if group_mean > grand_mean[marker]:
                n_clusters_high += 1
        f_statistic = (
            (between_sum_squares / (k - 1)) / (within_sum_squares / (n - k))
            if within_sum_squares > 0 and k > 1 and n > k
            else 0.0
        )
        results.append({"marker": marker, "F_statistic": f_statistic,
                        "n_clusters_high": n_clusters_high})
    return pd.DataFrame(results).sort_values("F_statistic", ascending=False)


def compute_cluster_marker_means(activity_df, markers, cluster_labels, method="highest_as"):
    """Compute mean marker activity in each cluster."""
    clusters = sorted(cluster_labels.unique())
    cluster_means = pd.DataFrame(index=clusters, columns=markers, dtype=float)
    for cluster in clusters:
        mask = cluster_labels == cluster
        cells = cluster_labels.index[mask]
        cluster_means.loc[cluster] = activity_df.loc[markers, cells].mean(axis=1).values
    if method == "highest_as":
        return cluster_means
    if method == "zscore":
        return cluster_means.apply(
            lambda row: (row - row.mean()) / (row.std() + 1e-10), axis=1
        )
    raise ValueError(f"Unknown method: {method}")


def assign_cell_types(cluster_means):
    """Assign each cluster to the marker with the highest mean activity."""
    return cluster_means.idxmax(axis=1), cluster_means.max(axis=1)
