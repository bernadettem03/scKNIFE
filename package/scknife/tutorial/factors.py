"""Consensus programs, representative nodes and plots for the fitted model."""

from types import SimpleNamespace
from ._plot_style import FIGURE_STYLE
from .core import output_dir
import json, hashlib, textwrap, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from matplotlib.patches import Patch


def row_z(x):
    x = np.log1p(np.asarray(x, dtype=float))
    return (x - x.mean(axis=1, keepdims=True)) / np.maximum(
        x.std(axis=1, keepdims=True), 1e-12
    )


def save(name):
    folder = output_dir() / "factor_interpretation_k4"
    folder.mkdir(exist_ok=True)
    plt.gcf().savefig(folder / (name + ".png"), dpi=170, bbox_inches="tight")


def heat(
    table, title, name, label="Row z-score of log(1 + mean activity)", signed=True
):
    fig, ax = plt.subplots(figsize=(6.3, max(3.0, 0.18 * len(table) + 1.45)))
    lim = max(float(np.max(np.abs(table.to_numpy()))), 1e-12)
    im = ax.imshow(
        table,
        aspect="auto",
        cmap="RdBu_r" if signed else "GnBu",
        vmin=-lim if signed else 0,
        vmax=lim,
    )
    ax.set_xticks(range(table.shape[1]), table.columns, rotation=40, ha="right")
    ax.set_yticks(range(len(table)), table.index)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=label)
    fig.tight_layout()
    save(name)
    return ax


def group_factors(a, graph, *, n_groups=4, bootstraps=200, seed=42):
    """Group fitted factors by bootstrapped W similarity and select a representative node dictionary."""
    OUT = output_dir() / "factor_interpretation_k4"
    OUT.mkdir(exist_ok=True)
    START = time.time()
    assert a.n_obs > 0 and a.obs_names.is_unique
    nodes = np.asarray(a.uns["node_names"]).astype(str)
    W = np.asarray(a.uns["W_graph"], dtype=float)
    hcols = sorted(
        [c for c in a.obs if c.startswith("H_graph_")],
        key=lambda c: int(c.rsplit("_", 1)[1]),
    )
    H = a.obs[hcols].to_numpy(float)
    N_FACTORS = H.shape[1]
    assert H.shape[0] == a.n_obs and W.shape == (len(nodes), N_FACTORS)
    assert (
        np.isfinite(W).all() and np.isfinite(H).all() and W.min() >= 0 and H.min() >= 0
    )
    meta = a.obs[["cellType", "patient_id", "timepoint"]].astype(str)
    types = sorted(meta.cellType.unique())
    G = graph
    kinds = np.array([str(G.nodes[n].get("kind", "")) for n in nodes])
    keep = (~np.isin(kinds, ["cell_type", "cell_state"])) & (W.sum(axis=1) > 0)
    W_ELIGIBLE = pd.DataFrame(
        W[keep], index=nodes[keep], columns=[str(i) for i in range(N_FACTORS)]
    )
    NODE_METADATA_ELIGIBLE = pd.DataFrame(
        {"node_name": nodes[keep], "node_kind": kinds[keep]}
    )
    factor_means = (
        pd.DataFrame(H, index=a.obs_names, columns=hcols)
        .groupby(meta.cellType)
        .mean()
        .reindex(types)
        .T
    )
    factor_means.to_csv(OUT / "factor_usage.csv")
    K_VALUES = range(2, min(13, N_FACTORS))
    W_Z = row_z(W_ELIGIBLE)
    rng = np.random.default_rng(seed)
    sums = {k: np.zeros((N_FACTORS, N_FACTORS)) for k in K_VALUES}
    for b in range(bootstraps):
        rows = rng.choice(len(W_Z), int(0.8 * len(W_Z)), replace=True)
        d = np.clip(1 - np.corrcoef(W_Z[rows].T), 0, 2)
        np.fill_diagonal(d, 0)
        for k in K_VALUES:
            lab = AgglomerativeClustering(
                n_clusters=k, metric="precomputed", linkage="average"
            ).fit_predict(d)
            sums[k] += lab[:, None] == lab[None, :]
    records = []
    solutions = {}
    for k in K_VALUES:
        c = sums[k] / bootstraps
        d = 1 - c
        np.fill_diagonal(d, 0)
        tree = linkage(squareform(d, checks=False), method="average")
        lab = fcluster(tree, t=k, criterion="maxclust") - 1
        sizes = np.bincount(lab)
        upper = c[np.triu_indices(N_FACTORS, 1)]
        records.append(
            {
                "k": k,
                "actual_groups": len(sizes),
                "PAC": np.mean((upper > 0.1) & (upper < 0.9)),
                "silhouette": silhouette_score(d, lab, metric="precomputed"),
                "min_size": sizes.min(),
                "max_size": sizes.max(),
            }
        )
        solutions[k] = (c, lab, tree)
    metrics = pd.DataFrame(records)
    # Review consensus diagnostics, then use the chosen four-group resolution.
    eligible = metrics[(metrics.min_size >= 2) & (metrics.actual_groups == metrics.k)]
    assert not eligible.empty
    assert n_groups in eligible.k.values
    metrics.to_csv(OUT / "k_diagnostics.csv", index=False)
    consensus, raw, tree = solutions[n_groups]
    raw_order = sorted(set(raw), key=lambda l: np.flatnonzero(raw == l)[0])
    labels = np.array([raw_order.index(l) for l in raw])
    groups = [f"CC{i + 1}" for i in range(n_groups)]
    ngroups = [f"NG{i + 1}" for i in range(n_groups)]
    COMPONENT_ASSIGNMENTS = pd.DataFrame(
        {"component": [str(i) for i in range(N_FACTORS)], "consensus_cluster": labels}
    )
    COMPONENT_ASSIGNMENTS.to_csv(OUT / "factor_groups.csv", index=False)
    order = np.lexsort((np.arange(N_FACTORS), labels))
    ticks = [f"H{i} · {groups[labels[i]]}" for i in order]
    DECAY_TAU = 5.0
    DECAY_CUMULATIVE_MASS = 0.95
    MAX_NODES_PER_CLUSTER_KIND = 20
    TOP_N = 20
    from scknife.tutorial_helpers import (
        score_nodes,
        assign_nodes_to_strongest_cluster,
        rank_nodes_with_exponential_decay,
        select_mixed_decay_budget,
    )

    ALL_NODE_SCORES = score_nodes(
        W_Z, NODE_METADATA_ELIGIBLE, COMPONENT_ASSIGNMENTS, n_groups
    )
    assigned = assign_nodes_to_strongest_cluster(ALL_NODE_SCORES)
    ranked, decay_audit = rank_nodes_with_exponential_decay(
        assigned, DECAY_TAU, DECAY_CUMULATIVE_MASS, MAX_NODES_PER_CLUSTER_KIND
    )
    local, local_audit = rank_nodes_with_exponential_decay(
        ALL_NODE_SCORES[ALL_NODE_SCORES.enrichment_delta_z > 0],
        DECAY_TAU,
        DECAY_CUMULATIVE_MASS,
        MAX_NODES_PER_CLUSTER_KIND,
    )
    selected = select_mixed_decay_budget(ranked, local, TOP_N, n_groups)
    assert selected.node_name.is_unique and len(selected) == TOP_N * n_groups
    assert selected.enrichment_delta_z.gt(0).all()
    selected.to_csv(OUT / "selected_nodes.csv", index=False)
    decay_audit.to_csv(OUT / "decay_audit.csv", index=False)
    wmean = np.column_stack(
        [W_ELIGIBLE.to_numpy()[:, labels == g].mean(axis=1) for g in range(n_groups)]
    )
    wz = pd.DataFrame(row_z(wmean), index=W_ELIGIBLE.index, columns=groups)
    profiles = np.vstack(
        [
            W[
                np.isin(
                    nodes, selected.loc[selected.consensus_cluster == g, "node_name"]
                )
            ].mean(axis=0)
            for g in range(n_groups)
        ]
    )
    node_raw = H @ profiles.T
    factor_raw = np.column_stack(
        [H[:, labels == g].mean(axis=1) for g in range(n_groups)]
    )
    node_mean = (
        pd.DataFrame(node_raw, index=a.obs_names, columns=ngroups)
        .groupby(meta.cellType)
        .mean()
        .reindex(types)
        .T
    )
    factor_mean = (
        pd.DataFrame(factor_raw, index=a.obs_names, columns=groups)
        .groupby(meta.cellType)
        .mean()
        .reindex(types)
        .T
    )
    node_type = pd.DataFrame(row_z(node_mean), index=ngroups, columns=types)
    factor_type = pd.DataFrame(row_z(factor_mean), index=groups, columns=types)
    node_mean.to_csv(OUT / "node_group_means.csv")
    factor_mean.to_csv(OUT / "factor_group_means.csv")
    summary = {
        "cells": a.n_obs,
        "factors": N_FACTORS,
        "selected_k": n_groups,
        "group_sizes": np.bincount(labels).tolist(),
        "bootstraps": bootstraps,
        "selected_nodes": len(selected),
        "seconds": time.time() - START,
        "cell_id_sha256": hashlib.sha256("\n".join(a.obs_names).encode()).hexdigest(),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(summary, flush=True)
    return SimpleNamespace(
        a=a,
        nodes=nodes,
        W=W,
        H=H,
        hcols=hcols,
        N_FACTORS=N_FACTORS,
        meta=meta,
        types=types,
        W_ELIGIBLE=W_ELIGIBLE,
        W_Z=W_Z,
        metrics=metrics,
        eligible=eligible,
        K_VALUES=K_VALUES,
        n_groups=n_groups,
        solutions=solutions,
        consensus=consensus,
        labels=labels,
        groups=groups,
        ngroups=ngroups,
        order=order,
        ticks=ticks,
        selected=selected,
        ranked=ranked,
        entity_kinds=sorted(ranked.loc[ranked.decay_selected, "node_kind"].unique()),
        wz=wz,
        node_raw=node_raw,
        node_type=node_type,
        factor_type=factor_type,
        factor_means=factor_means,
    )


@plt.rc_context(FIGURE_STYLE)
def plot_usage(programs):
    """Show mean H usage by source cell type."""
    factor_means = programs.factor_means
    return heat(
        factor_means,
        "Mean H activity by source cell type",
        "factor_usage",
        "Mean H activity",
        False,
    )


@plt.rc_context(FIGURE_STYLE)
def plot_group_selection(programs):
    """Compare ambiguity and silhouette across candidate group counts."""
    K_VALUES = programs.K_VALUES
    eligible = programs.eligible
    metrics = programs.metrics
    n_groups = programs.n_groups
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.15))
    for ax, col, title in zip(
        axes,
        ["PAC", "silhouette"],
        [
            "Ambiguous fraction (lower is better)",
            "Silhouette (higher is better)",
        ],
    ):
        ax.plot(metrics.k, metrics[col], "o-", color="#247b87")
        ax.scatter(
            eligible.k,
            eligible[col],
            s=110,
            facecolors="none",
            edgecolors="#d7832c",
            label="No singleton groups",
        )
        ax.axvline(n_groups, color="#9b456b", ls="--", label=f"Selected k={n_groups}")
        ax.set_xlabel("Number of factor groups (k)")
        ax.set_title(title, fontsize=10)
        ax.set_xticks(list(K_VALUES))
    axes[0].legend(fontsize=9)
    fig.tight_layout()
    save("choose_k")
    return axes


@plt.rc_context(FIGURE_STYLE)
def plot_consensus(programs):
    """Show co-clustering frequency across bootstrap runs."""
    consensus = programs.consensus
    order = programs.order
    ticks = programs.ticks
    return heat(
        pd.DataFrame(consensus[np.ix_(order, order)], index=ticks, columns=ticks),
        "Bootstrap factor co-clustering",
        "consensus",
        "Co-clustering fraction",
        False,
    )


@plt.rc_context(FIGURE_STYLE)
def plot_resolutions(programs):
    """Compare membership at two, three and four groups."""
    N_FACTORS = programs.N_FACTORS
    solutions = programs.solutions
    sensitivity = []
    for k in [2, 3, 4]:
        lab = solutions[k][1]
        ordered_labels = sorted(set(lab), key=lambda v: np.flatnonzero(lab == v)[0])
        sensitivity.append([ordered_labels.index(v) + 1 for v in lab])
    fig, ax = plt.subplots(figsize=(7.0, 2.2))
    ax.imshow(sensitivity, cmap="tab10", aspect="auto", vmin=1, vmax=10)
    for i, row in enumerate(sensitivity):
        for j, v in enumerate(row):
            ax.text(j, i, str(v), ha="center", va="center", color="white")
    ax.set_xticks(range(N_FACTORS), [f"H{i}" for i in range(N_FACTORS)], rotation=45)
    ax.set_yticks(range(3), ["k=2", "k=3", "k=4"])
    ax.set_title("Factor membership at nearby resolutions")
    fig.tight_layout()
    save("resolution_sensitivity")
    return ax


@plt.rc_context(FIGURE_STYLE)
def plot_node_dictionary(programs, by="factor"):
    """Show all selected nodes across factor groups or individual factors."""
    W_ELIGIBLE = programs.W_ELIGIBLE
    W_Z = programs.W_Z
    groups = programs.groups
    labels = programs.labels
    n_groups = programs.n_groups
    ngroups = programs.ngroups
    order = programs.order
    selected = programs.selected
    wz = programs.wz
    rows = selected.sort_values(["consensus_cluster", "cluster_rank"])
    row_groups = rows.consensus_cluster.to_numpy()
    row_names = [
        textwrap.fill(f"{ngroups[g]}: {n} [{kind}]", 52)
        for g, n, kind in zip(row_groups, rows.node_name, rows.node_kind)
    ]

    def combined_node_heatmap(
        values, columns, title, filename, color_label, factor_columns=False
    ):
        fig, ax = plt.subplots(figsize=(8.8, max(4.2, 0.30 * len(rows) + 1.5)))
        limit = max(float(np.abs(values).max()), 1e-12)
        im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
        ax.set_yticks(range(len(rows)), row_names, fontsize=9)
        ax.set_xticks(range(len(columns)), columns, rotation=45, ha="right")
        for boundary in np.flatnonzero(np.diff(row_groups) != 0) + 0.5:
            ax.axhline(boundary, color="#273e50", lw=1.2)
        if factor_columns:
            for boundary in np.flatnonzero(np.diff(labels[order]) != 0) + 0.5:
                ax.axvline(boundary, color="#273e50", lw=1.2)
            top = ax.secondary_xaxis("top")
            top.set_xticks(
                [np.flatnonzero(labels[order] == g).mean() for g in range(n_groups)],
                groups,
            )
            top.tick_params(length=0)
        ax.set_title(title, pad=30 if factor_columns else 12)
        fig.colorbar(im, ax=ax, label=color_label, fraction=0.025, pad=0.03)
        fig.tight_layout()
        save(filename)
        return ax

    if by == "group":
        return combined_node_heatmap(
            wz.loc[rows.node_name].to_numpy(),
            groups,
            "Representative biological nodes · NG1–NG4",
            "node_dictionary_all",
            "Node-wise z-score of log(1 + mean W)",
        )
    if by != "factor":
        raise ValueError("Choose factor or group")
    return combined_node_heatmap(
        W_Z[W_ELIGIBLE.index.get_indexer(rows.node_name)][:, order],
        [f"H{i}" for i in order],
        "Biological nodes across factors · NG1–NG4",
        "node_factors_all",
        "Node-wise z-score of log(1 + W)",
        factor_columns=True,
    )


@plt.rc_context(FIGURE_STYLE)
def plot_composition(programs):
    """Count selected nodes by entity type and program."""
    ngroups = programs.ngroups
    selected = programs.selected
    composition = pd.crosstab(selected.consensus_cluster, selected.node_kind)
    composition.index = ngroups
    ax = composition.plot.bar(stacked=True, figsize=(6.3, 3.5), colormap="tab10", rot=0)
    ax.set_ylabel("Selected nodes")
    ax.set_xlabel("Node group")
    ax.set_title("Biological node-kind composition")
    ax.legend(
        title="Node kind", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=10
    )
    plt.tight_layout()
    save("composition")
    return ax


@plt.rc_context(FIGURE_STYLE)
def plot_entities(programs, kind="geneset"):
    """Show a single entity category across all programs."""
    N_FACTORS = programs.N_FACTORS
    W_ELIGIBLE = programs.W_ELIGIBLE
    W_Z = programs.W_Z
    groups = programs.groups
    labels = programs.labels
    n_groups = programs.n_groups
    ngroups = programs.ngroups
    order = programs.order
    ranked = programs.ranked
    rows = ranked[ranked.decay_selected & ranked.node_kind.eq(kind)]
    if rows.empty:
        raise ValueError(f"No selected nodes for {kind}")
    rows = rows.sort_values(["consensus_cluster", "category_rank"])
    values = W_Z[W_ELIGIBLE.index.get_indexer(rows.node_name)][:, order]
    row_groups = rows.consensus_cluster.to_numpy()
    names = [
        textwrap.fill(f"{ngroups[g]}: {n}", 52)
        for g, n in zip(row_groups, rows.node_name)
    ]
    fig, ax = plt.subplots(figsize=(8.8, max(4.2, 0.30 * len(rows) + 1.5)))
    limit = max(float(np.abs(values).max()), 1e-12)
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-limit, vmax=limit)
    ax.set_yticks(range(len(rows)), names, fontsize=9)
    ax.set_xticks(range(N_FACTORS), [f"H{i}" for i in order], rotation=45, ha="right")
    for boundary in np.flatnonzero(np.diff(row_groups) != 0) + 0.5:
        ax.axhline(boundary, color="#273e50", lw=1.2)
    for boundary in np.flatnonzero(np.diff(labels[order]) != 0) + 0.5:
        ax.axvline(boundary, color="#273e50", lw=1.2)
    top = ax.secondary_xaxis("top")
    centers = [np.flatnonzero(labels[order] == g).mean() for g in range(n_groups)]
    top.set_xticks(centers, groups)
    top.tick_params(length=0)
    ax.set_title(f"{kind.replace('_', ' ')} · all programs", pad=30)
    ax.set_xlabel("Factors ordered by factor group")
    fig.colorbar(
        im, ax=ax, label="Node-wise z-score of log(1 + W)", fraction=0.025, pad=0.03
    )
    fig.tight_layout()
    save(f"category_all_{kind}")
    return ax


@plt.rc_context(FIGURE_STYLE)
def plot_group_activity(programs, group="NG"):
    """Compare reconstructed node-group or factor-group activity."""
    factor_type = programs.factor_type
    node_type = programs.node_type
    if group == "NG":
        return heat(node_type, "Reconstructed node-group activity", "node_activity")
    if group == "CC":
        return heat(factor_type, "Shared W–H factor-group activity", "factor_activity")
    raise ValueError("Choose NG or CC")


@plt.rc_context(FIGURE_STYLE)
def plot_associations(programs):
    """Connect relative W and H profiles to cell types."""
    factor_type = programs.factor_type
    groups = programs.groups
    n_groups = programs.n_groups
    ngroups = programs.ngroups
    selected = programs.selected
    types = programs.types
    wz = programs.wz
    OUT = output_dir() / "factor_interpretation_k4"
    node_factor = pd.DataFrame(
        np.vstack(
            [
                wz.loc[selected.loc[selected.consensus_cluster == g, "node_name"]].mean(
                    axis=0
                )
                for g in range(n_groups)
            ]
        ),
        index=ngroups,
        columns=groups,
    )
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ys = [np.linspace(0, 1, len(x)) for x in [ngroups, groups, types]]
    edges = []
    for layer, table in enumerate([node_factor, factor_type]):
        for i, src in enumerate(table.index):
            for j, dst in enumerate(table.columns):
                v = float(table.iloc[i, j])
                edges.append({"source": src, "target": dst, "relative_score": v})
                ax.plot(
                    [layer, layer + 1],
                    [ys[layer][i], ys[layer + 1][j]],
                    color="#b24b55" if v >= 0 else "#357eaa",
                    lw=0.3 + 2 * abs(v),
                    alpha=0.38,
                )
    for x, names in enumerate([ngroups, groups, types]):
        ax.scatter([x] * len(names), ys[x], s=85, color="#273e50", zorder=3)
        for y, name in zip(ys[x], names):
            ax.text(x + 0.04, y, name, va="center", fontsize=9.5)
    for x, title in enumerate(
        ["W node groups", "Shared factor groups", "Source cell types"]
    ):
        ax.text(x, 1.08, title, weight="bold")
    ax.set_xlim(-0.05, 2.5)
    ax.set_ylim(-0.06, 1.15)
    ax.axis("off")
    ax.legend(
        handles=[
            Patch(color="#b24b55", label="Positive relative score"),
            Patch(color="#357eaa", label="Negative relative score"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=2,
    )
    plt.tight_layout()
    save("associations")
    pd.DataFrame(edges).to_csv(OUT / "all_associations.csv", index=False)
    return ax


@plt.rc_context(FIGURE_STYLE)
def plot_cell_profiles(programs):
    """Show node-group activity across ordered cells with metadata tracks."""
    a = programs.a
    meta = programs.meta
    n_groups = programs.n_groups
    ngroups = programs.ngroups
    node_raw = programs.node_raw
    logged = np.log1p(node_raw)
    cell_z = (logged - logged.mean(axis=0)) / np.maximum(logged.std(axis=0), 1e-12)
    ordered = meta.sort_values(
        ["cellType", "patient_id", "timepoint"], kind="stable"
    ).index
    positions = a.obs_names.get_indexer(ordered)
    fig, axes = plt.subplots(
        4,
        1,
        figsize=(7.1, 4.7),
        gridspec_kw={"height_ratios": [0.25, 0.25, 0.25, 3]},
        layout="constrained",
    )
    legend = []
    for ax, key in zip(axes[:3], ["cellType", "patient_id", "timepoint"]):
        cats = sorted(meta[key].unique())
        codes = pd.Categorical(meta.loc[ordered, key], categories=cats).codes
        ax.imshow(codes[None, :], aspect="auto", cmap="tab20", vmin=0, vmax=19)
        ax.set_yticks([0], [key])
        ax.set_xticks([])
        for j, c in enumerate(cats):
            legend.append(Patch(color=plt.get_cmap("tab20")(j), label=f"{key}: {c}"))
    im = axes[-1].imshow(
        cell_z[positions].T, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5
    )
    axes[-1].set_yticks(range(n_groups), ngroups)
    axes[-1].set_xlabel(
        f"All {a.n_obs:,} cells · ordered by cell type, patient and timepoint"
    )
    fig.colorbar(im, ax=axes[-1], label="Per-group z-score of log(1 + activity)")
    fig.legend(handles=legend, loc="outside lower center", ncol=4, fontsize=9)
    save("all_cells")
    return axes


@plt.rc_context(FIGURE_STYLE)
def plot_coverage(programs):
    """Show cell counts per patient and cell type."""
    meta = programs.meta
    types = programs.types
    coverage = pd.crosstab(meta.cellType, meta.patient_id).reindex(types)
    return heat(
        coverage,
        "Cells available per patient and cell type",
        "coverage",
        "Number of cells",
        False,
    )


@plt.rc_context(FIGURE_STYLE)
def plot_constituents(programs):
    """Show three ranked constituents of the first node group."""
    H = programs.H
    W = programs.W
    meta = programs.meta
    nodes = programs.nodes
    selected = programs.selected
    types = programs.types
    chosen = selected[selected.consensus_cluster == 0].head(3).node_name.tolist()
    fig, axes = plt.subplots(1, 3, figsize=(8.1, 4.0))
    for ax, node in zip(axes, chosen):
        score = np.log1p(H @ W[np.flatnonzero(nodes == node)[0]])
        ax.violinplot(
            [score[meta.cellType.to_numpy() == t] for t in types], showmedians=True
        )
        ax.set_xticks(range(1, len(types) + 1), types, rotation=60, ha="right")
        ax.set_title(textwrap.fill(node, 28), fontsize=9.5)
        ax.set_ylabel("log(1 + reconstructed node activity)")
    plt.tight_layout()
    save("constituents")
    return axes
