"""Node-selection helpers for the factor-interpretation tutorial."""
import numpy as np
import pandas as pd

def score_nodes(node_zscores, metadata, assignments, selected_k):
    labels_local = assignments["consensus_cluster"].to_numpy()
    records = []
    for cluster_id in range(selected_k):
        inside = labels_local == cluster_id
        outside = ~inside
        frame = metadata.copy()
        frame["consensus_cluster"] = cluster_id
        frame["n_cluster_components"] = int(inside.sum())
        frame["mean_z_inside"] = node_zscores[:, inside].mean(axis=1)
        frame["mean_z_outside"] = node_zscores[:, outside].mean(axis=1)
        frame["enrichment_delta_z"] = (
            frame["mean_z_inside"] - frame["mean_z_outside"]
        )
        frame["category_percentile"] = frame.groupby(
            "node_kind"
        )["enrichment_delta_z"].rank(method="average", pct=True)
        records.append(frame)
    return pd.concat(records, ignore_index=True)

def assign_nodes_to_strongest_cluster(scores):
    """Assign each node once to its strongest positive-enrichment cluster."""
    ordered = scores.sort_values(
        [
            "node_name",
            "enrichment_delta_z",
            "mean_z_inside",
            "consensus_cluster",
        ],
        ascending=[True, False, False, True],
        kind="mergesort",
    )
    assigned = ordered.drop_duplicates("node_name", keep="first").copy()
    assigned = assigned[assigned["enrichment_delta_z"] > 0].copy()
    return assigned.sort_values(
        ["consensus_cluster", "node_kind", "node_name"],
        kind="mergesort",
    ).reset_index(drop=True)

def rank_nodes_with_exponential_decay(
    assigned,
    tau=5.0,
    cumulative_mass=0.95,
    max_per_group=20,
):
    """Rank independently within cluster × node kind and select a decay prefix."""
    if tau <= 0:
        raise ValueError("DECAY_TAU must be positive")
    if not 0 < cumulative_mass <= 1:
        raise ValueError("DECAY_CUMULATIVE_MASS must be in (0, 1]")
    if int(max_per_group) < 1:
        raise ValueError("MAX_NODES_PER_CLUSTER_KIND must be positive")

    ranked_parts = []
    summary_records = []
    group_columns = ["consensus_cluster", "node_kind"]
    for (cluster_id, node_kind), group in assigned.groupby(
        group_columns, sort=True, observed=True
    ):
        group = group.sort_values(
            ["enrichment_delta_z", "mean_z_inside", "node_name"],
            ascending=[False, False, True],
            kind="mergesort",
        ).reset_index(drop=True)
        group["category_rank"] = np.arange(1, len(group) + 1)
        group["decay_weight"] = np.exp(
            -(group["category_rank"].to_numpy(dtype=float) - 1.0) / tau
        )
        total_weight = float(group["decay_weight"].sum())
        group["decay_mass_fraction"] = group["decay_weight"] / total_weight
        group["decay_cumulative_mass"] = group[
            "decay_mass_fraction"
        ].cumsum()

        reached = np.flatnonzero(
            group["decay_cumulative_mass"].to_numpy() >= cumulative_mass
        )
        uncapped_cutoff = int(reached[0] + 1) if len(reached) else len(group)
        cutoff_rank = min(uncapped_cutoff, int(max_per_group))
        group["decay_cutoff_rank"] = cutoff_rank
        group["category_quota"] = cutoff_rank
        group["decay_selected"] = group["category_rank"] <= cutoff_rank
        group["decay_backfill"] = False
        group["selection_reason"] = np.where(
            group["decay_selected"], "decay_prefix", "not_selected"
        )
        achieved_mass = float(
            group.loc[cutoff_rank - 1, "decay_cumulative_mass"]
        )
        summary_records.append({
            "consensus_cluster": int(cluster_id),
            "node_kind": str(node_kind),
            "eligible_positive": int(len(group)),
            "decay_selected": int(cutoff_rank),
            "decay_tau": float(tau),
            "target_cumulative_mass": float(cumulative_mass),
            "achieved_cumulative_mass": achieved_mass,
            "maximum_per_cluster_kind": int(max_per_group),
            "cap_reached": bool(uncapped_cutoff > max_per_group),
            "target_reached": bool(achieved_mass >= cumulative_mass),
        })
        ranked_parts.append(group)

    ranked = pd.concat(ranked_parts, ignore_index=True)
    ranked = ranked.sort_values(
        ["consensus_cluster", "node_kind", "category_rank"],
        kind="mergesort",
    ).reset_index(drop=True)
    summary = pd.DataFrame(summary_records).sort_values(
        ["consensus_cluster", "node_kind"], kind="mergesort"
    ).reset_index(drop=True)
    return ranked, summary

def select_mixed_decay_budget(
    assigned_ranked, cluster_local_ranked, per_cluster_budget, selected_k
):
    """Fill a globally unique fixed budget, prioritizing strongest assignments."""
    selected_parts = []
    used_nodes = set()
    priority_columns = [
        "decay_weight",
        "enrichment_delta_z",
        "mean_z_inside",
        "node_name",
    ]
    priority_ascending = [False, False, False, True]

    assigned_prefix_counts = (
        assigned_ranked[assigned_ranked["decay_selected"]]
        .groupby("consensus_cluster")
        .size()
        .reindex(range(selected_k), fill_value=0)
    )
    cluster_order = assigned_prefix_counts.sort_values(
        ascending=True, kind="mergesort"
    ).index.tolist()

    for cluster_id in cluster_order:
        assigned_rows = assigned_ranked[
            assigned_ranked["consensus_cluster"].eq(cluster_id)
            & assigned_ranked["decay_selected"]
            & ~assigned_ranked["node_name"].isin(used_nodes)
        ].sort_values(
            priority_columns,
            ascending=priority_ascending,
            kind="mergesort",
        )
        chosen = assigned_rows.head(per_cluster_budget).copy()
        chosen["decay_backfill"] = False
        chosen["selection_reason"] = "strongest_cluster_decay_prefix"

        if len(chosen) < per_cluster_budget:
            local_rows = cluster_local_ranked[
                cluster_local_ranked["consensus_cluster"].eq(cluster_id)
                & cluster_local_ranked["decay_selected"]
                & ~cluster_local_ranked["node_name"].isin(used_nodes)
                & ~cluster_local_ranked["node_name"].isin(
                    chosen["node_name"]
                )
            ].sort_values(
                priority_columns,
                ascending=priority_ascending,
                kind="mergesort",
            )
            local_fill = local_rows.head(
                per_cluster_budget - len(chosen)
            ).copy()
            local_fill["decay_backfill"] = True
            local_fill["selection_reason"] = "cluster_local_decay_prefix"
            chosen = pd.concat([chosen, local_fill], ignore_index=True)

        if len(chosen) < per_cluster_budget:
            remaining = cluster_local_ranked[
                cluster_local_ranked["consensus_cluster"].eq(cluster_id)
                & ~cluster_local_ranked["node_name"].isin(used_nodes)
                & ~cluster_local_ranked["node_name"].isin(
                    chosen["node_name"]
                )
            ].sort_values(
                priority_columns,
                ascending=priority_ascending,
                kind="mergesort",
            )
            backfill = remaining.head(
                per_cluster_budget - len(chosen)
            ).copy()
            backfill["decay_backfill"] = True
            backfill["selection_reason"] = "budget_backfill"
            chosen = pd.concat([chosen, backfill], ignore_index=True)

        if len(chosen) != per_cluster_budget:
            raise RuntimeError(
                f"Cluster {cluster_id} has only {len(chosen)} available "
                f"unique positive nodes for a {per_cluster_budget}-node budget"
            )
        chosen = chosen.sort_values(
            priority_columns,
            ascending=priority_ascending,
            kind="mergesort",
        ).reset_index(drop=True)
        chosen["cluster_rank"] = np.arange(1, len(chosen) + 1)
        used_nodes.update(chosen["node_name"].tolist())
        selected_parts.append(chosen)

    return pd.concat(selected_parts, ignore_index=True).sort_values(
        ["consensus_cluster", "cluster_rank"], kind="mergesort"
    ).reset_index(drop=True)
