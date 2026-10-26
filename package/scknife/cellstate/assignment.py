"""Activity-based cell-state assignment helpers."""

import pandas as pd


def compute_cluster_state_means(adata, found_states, cluster_col="leiden_clusters"):
    """Assign each cluster to its highest-activity state marker."""
    cluster_means = pd.DataFrame()
    state_assignment = {}
    for cluster in sorted(adata.obs[cluster_col].unique()):
        mask = adata.obs[cluster_col] == cluster
        state_means = adata.obs.loc[mask, found_states].mean()
        state_means_sorted = state_means.sort_values(ascending=False)
        state_assignment[cluster] = {
            "cell_state": state_means_sorted.index[0],
            "score": state_means_sorted.iloc[0],
            "n_cells": int(mask.sum()),
            "all_scores": state_means_sorted,
        }
        cluster_means[f"Cluster_{cluster}"] = state_means
    return state_assignment, cluster_means


def build_state_mapping_table(state_assignment):
    """Convert a cluster-state assignment dictionary to a tidy table."""
    return pd.DataFrame([
        {"Cluster": cluster,
         "Assigned_CellState": state_assignment[cluster]["cell_state"],
         "Mean_Score": state_assignment[cluster]["score"],
         "N_Cells": state_assignment[cluster]["n_cells"]}
        for cluster in sorted(state_assignment)
    ])
