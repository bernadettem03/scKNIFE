"""Activity-weighted sender–receiver signaling on a saved fit."""

from contextlib import contextmanager, ExitStack
from pathlib import Path
from types import SimpleNamespace
import os, hashlib, json, textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from statsmodels.stats.multitest import multipletests
from scknife.signaling import activity_lr as lr, common, plotting as plots
from scknife.tutorial_inputs import prepare_signaling_input
from ._plot_style import FIGURE_STYLE
from .core import load_fit, output_dir, figure, repository_root


# Match the cell-type colors across categorical views of the same network.
_CELL_COLORS = {
    "B_cell": "#569bc0", "Cancer_cell": "#d98565", "Endothelial_cell": "#62a69d",
    "Fibroblast": "#b49a57", "Mast_cell": "#b887a9", "Myeloid_cell": "#8190b5",
    "T_cell": "#82a76d", "pDC": "#a18cbb",
}



@contextmanager
def _settings(module, **values):
    """Apply legacy module settings only for the enclosed call."""
    previous = {key: getattr(module, key) for key in values}
    try:
        for key, value in values.items():
            setattr(module, key, value)
        yield
    finally:
        for key, value in previous.items():
            setattr(module, key, value)


def score_signaling(fit=None, *, activity_path=None, groupby="cellType", interaction_table=None):
    """Return pooled scores and interaction details grouped by the selected metadata column."""
    base = output_dir()
    out = base / "signaling"
    out.mkdir(exist_ok=True)
    adata = load_fit() if fit is None else fit
    activity_path = base / "activity_scores.parquet" if activity_path is None else Path(activity_path)
    metadata = adata.obs
    key = groupby
    db_path = prepare_signaling_input(
        repository_root() / "tutorial/inputs",
        interaction_table or os.environ.get("SCKNIFE_INTERACTION_DB"),
    )
    with ExitStack() as stack:
        for module in (lr, common):
            stack.enter_context(
                _settings(
                    module,
                    ACTIVITY_PATH=activity_path,
                    CELLCHAT_DB_PATH=db_path,
                )
            )
        assert key in metadata and not metadata[key].isna().any()
        db, wanted = lr.load_interactions()
        _, nodes, cells = lr.activity_schema()
        if cells != list(adata.obs_names):
            raise ValueError("Activity columns must match the fitted cells in order.")
        values, genes, coverage = lr.stream_activity_subset(wanted, nodes, cells, 256)
        groups = sorted(metadata[key].astype(str).unique())
        assert len(groups) > 1
        labels = (
            metadata[key]
            .astype(str)
            .map({g: i for i, g in enumerate(groups)})
            .to_numpy(dtype=np.int16)
        )
        means, counts = lr.group_means(values, labels, len(groups))
        specificity = means / (means.mean(axis=1)[:, None] + 1e-10)
        pairs = lr.build_pairs(db, genes, means)
        lr.validate_vectorization(means, specificity, pairs)
        comm = lr.communication_matrix(means, specificity, pairs)
        assert np.isfinite(comm).all() and (comm >= 0).all()
        details = lr.interaction_details(means, specificity, pairs, groups)
        reconstructed = details.pivot_table(
            index="source",
            columns="target",
            values="weighted_score",
            aggfunc="sum",
            fill_value=0,
        ).reindex(index=groups, columns=groups, fill_value=0)
        assert np.allclose(comm, reconstructed.to_numpy(), rtol=1e-10, atol=1e-10)
        pd.DataFrame(comm, index=groups, columns=groups).to_csv(
            out / "communication_matrix.csv"
        )
        details.to_csv(
            out / "interaction_details.csv.gz", index=False, compression="gzip"
        )
        pd.DataFrame({"reference_cell_type": groups, "cells": counts}).to_csv(
            out / "group_counts.csv", index=False
        )
        pairs.to_csv(out / "matched_interaction_components.csv", index=False)
    return SimpleNamespace(
        groups=groups,
        labels=labels,
        values=values,
        specificity=specificity,
        pairs=pairs,
        matrix=pd.DataFrame(comm, index=groups, columns=groups),
        details=details,
        out=out,
    )


@plt.rc_context(FIGURE_STYLE)
def plot_matrix(signals):
    """Plot the total signaling score for every sender–receiver pair."""
    comm, groups, out = signals.matrix, signals.groups, signals.out
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    im = ax.imshow(comm, cmap="GnBu", aspect="equal")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([textwrap.fill(g, 16) for g in groups], rotation=35, ha="right")
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([textwrap.fill(g, 20) for g in groups])
    ax.set_xlabel("Receiver cell type")
    ax.set_ylabel("Sender cell type")
    ax.set_title("Sender–receiver activity")
    fig.colorbar(im, ax=ax, label="Weighted activity score")
    fig.tight_layout()
    fig.savefig(out / "signaling_heatmap.png", dpi=160, bbox_inches="tight")
    return ax


@plt.rc_context(FIGURE_STYLE)
def plot_contributions(signals):
    """Show the strongest ligand–receptor components across selected group pairs."""
    details, out = signals.details, signals.out
    totals = (
        details.groupby(["ligand", "receptor"], sort=True)
        .weighted_score.sum()
        .sort_values(ascending=False, kind="stable")
        .head(10)
    )
    contexts = (
        details.groupby(["source", "target"], sort=True)
        .weighted_score.sum()
        .sort_values(ascending=False, kind="stable")
        .head(6)
    )
    agg = details.groupby(
        ["ligand", "receptor", "source", "target"], sort=True
    ).weighted_score.sum()
    rows = []
    for y, (lig, rec) in enumerate(totals.index):
        for x, (sender, receiver) in enumerate(contexts.index):
            score = float(agg.get((lig, rec, sender, receiver), 0))
            rows.append(
                dict(
                    ligand=lig,
                    receptor=rec,
                    source=sender,
                    target=receiver,
                    score=score,
                    x=x,
                    y=y,
                )
            )
    bubbles = pd.DataFrame(rows)
    bubbles.to_csv(out / "selected_lr_bubbles.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    positive = bubbles[bubbles.score > 0]
    sizes = 185 * positive.score / positive.score.max()
    scat = ax.scatter(
        positive.x,
        positive.y,
        s=sizes,
        c=np.log1p(positive.score),
        cmap="GnBu",
        edgecolors="#526777",
        linewidths=0.3,
    )
    ax.set_xticks(range(len(contexts)))
    ax.set_xticklabels(
        [
            a.replace("_cell", " cell") + " → " + b.replace("_cell", " cell")
            for a, b in contexts.index
        ],
        rotation=35, ha="right", fontsize=8.5,
    )
    ax.set_yticks(range(len(totals)))
    ax.set_yticklabels([a + " → " + b for a, b in totals.index], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("Sender → receiver")
    ax.set_ylabel("Ligand → receptor component")
    ax.set_title("Ligand–receptor contributions")
    fig.colorbar(scat, ax=ax, label="log(1 + weighted score)")
    fig.tight_layout()
    fig.savefig(out / "signaling_lr_bubbles.png", dpi=160, bbox_inches="tight")
    return ax


def permutation_fdr(signals, n_permutations=1000):
    """Return BH-adjusted p-values for global sender–receiver totals."""
    out = signals.out / "inference"
    out.mkdir(exist_ok=True)
    labels, groups = signals.labels, signals.groups
    sig = hashlib.sha256(labels.tobytes() + json.dumps(groups).encode()).hexdigest()
    with _settings(lr, DATA_DIR=out):
        pvalues, completed = lr.run_permutations(
            signals.values,
            labels,
            signals.matrix.to_numpy(),
            signals.specificity,
            signals.pairs,
            len(groups),
            n_permutations,
            50,
            False,
            sig,
        )
    fdr = pd.DataFrame(
        multipletests(pvalues.ravel(), method="fdr_bh")[1].reshape(pvalues.shape),
        index=groups,
        columns=groups,
    )
    pd.DataFrame(pvalues, index=groups, columns=groups).to_csv(
        out / "communication_pvalues.csv"
    )
    fdr.to_csv(out / "communication_fdr.csv")
    fdr.attrs["permutations"] = completed
    return fdr


@plt.rc_context(FIGURE_STYLE)
def plot_fdr(fdr):
    """Plot adjusted permutation significance by sender and receiver."""
    groups = fdr.index
    out = output_dir() / "signaling/inference"

    def heatmap(v, title, filename, label, diverging=False):
        fig, ax = plt.subplots(figsize=(5.6, 4.6))
        cmap = plt.get_cmap("RdBu_r" if diverging else "Purples").copy()
        cmap.set_bad("#dddddd")
        bound = max(float(np.nanmax(np.abs(v))), 1e-12)
        im = ax.imshow(
            v, cmap=cmap, aspect="equal", vmin=-bound if diverging else 0, vmax=bound
        )
        ax.set_xticks(range(len(groups)), groups, rotation=40, ha="right")
        ax.set_yticks(range(len(groups)), groups)
        ax.set_xlabel("Receiver cell type")
        ax.set_ylabel("Sender cell type")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label=label)
        fig.tight_layout()
        fig.savefig(out / filename, dpi=170)
        return ax

    return heatmap(
        -np.log10(fdr),
        f"Cell-label permutation FDR · {fdr.attrs['permutations']:,} permutations",
        "permutation_fdr.png",
        "−log10(BH-adjusted p-value)",
    )


def centrality(signals):
    """Return weighted network centrality after excluding self-links."""
    observed, groups = signals.matrix.to_numpy(), signals.groups
    out = signals.out / "inference"
    out.mkdir(exist_ok=True)
    g = nx.DiGraph()
    g.add_nodes_from(groups)
    for i, a in enumerate(groups):
        for j, b in enumerate(groups):
            if i != j and observed[i, j] > 0:
                g.add_edge(
                    a,
                    b,
                    weight=float(observed[i, j]),
                    distance=1 / float(observed[i, j]),
                )
    between = nx.betweenness_centrality(g, weight="distance", normalized=True)
    rank = nx.pagerank(g, weight="weight")
    central = pd.DataFrame(
        {
            "cell_type": groups,
            "outgoing_strength": [g.out_degree(a, weight="weight") for a in groups],
            "incoming_strength": [g.in_degree(a, weight="weight") for a in groups],
            "betweenness": [between[a] for a in groups],
            "pagerank": [rank[a] for a in groups],
        }
    ).set_index("cell_type")
    central.to_csv(out / "network_centrality.csv")
    return central


@plt.rc_context(FIGURE_STYLE)
def plot_centrality(central):
    """Show incoming and outgoing strength, betweenness and PageRank."""
    out = output_dir() / "signaling/inference"
    groups = central.index
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0), sharey=True)
    with _settings(plots, H_COLORS={**plots.H_COLORS, **_CELL_COLORS}):
        colors = plots.group_colors(list(groups))
    for ax, col in zip(axes.flat, central):
        ax.barh(groups, central[col], height=0.62, color=colors)
        ax.set_title(col.replace("_", " ").title())
        ax.set_xlabel("Score")
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    axes[0, 0].invert_yaxis()
    fig.suptitle("Network centrality")
    fig.tight_layout()
    fig.savefig(out / "network_centrality.png", dpi=170)
    return axes


def pathway_scores(signals, pathway="VEGF"):
    """Return pathway-specific sender–receiver totals and ligand–receptor detail."""
    SELECTED_PATHWAY = pathway
    details, groups = signals.details, signals.groups
    frame = details.loc[details.pathway_name.astype(str) == SELECTED_PATHWAY].copy()
    if frame.empty:
        raise ValueError("No interactions for " + SELECTED_PATHWAY)
    path_comm = frame.pivot_table(
        index="source",
        columns="target",
        values="weighted_score",
        aggfunc="sum",
        fill_value=0,
    ).reindex(index=groups, columns=groups, fill_value=0)
    return path_comm, frame


@plt.rc_context(FIGURE_STYLE)
def _render(fn, args, folder, filename, **kwargs):
    dest = output_dir() / folder
    with _settings(
        plots, FIGSIZE=(5.6, 4.6), WIDE_FIGSIZE=(7.0, 3.5), NODE_LABEL_SIZE_PT=8,
        H_COLORS={**plots.H_COLORS, **_CELL_COLORS}, COLORMAP_SEQUENTIAL="GnBu",
        EDGE_ALPHA=0.45, TOP_LR_PAIRS=15, OUTPUT_DPI=180,
    ):
        fn(*args, dest, **kwargs)
    return figure(folder / filename, width=600 if fn is plots.plot_sender_receiver else 490)


def plot_strength(signals):
    """Display outgoing and incoming scores, including self-links."""
    return _render(
        plots.plot_sender_receiver,
        (signals.matrix,),
        Path("signaling/general"),
        "01_sender_receiver.png",
    )


def plot_network(signals):
    """Display the strongest directional links in the pooled network."""
    return _render(
        plots.plot_network,
        (signals.matrix,),
        Path("signaling/general"),
        "06_hierarchical_network.png",
        title="Cell-type communication network",
    )


def plot_pathway(matrix, pathway="VEGF"):
    """Display a heatmap with a pathway-specific color scale."""
    return _render(
        plots.plot_heatmaps,
        (matrix, None),
        Path("signaling/pathways") / pathway,
        "03_communication_heatmap.png",
        title=f"{pathway} communication",
    )


def plot_pathway_pairs(details, pathway="VEGF"):
    """Rank up to 15 ligand–receptor pairs in a chosen pathway."""
    return _render(
        plots.plot_top_lr,
        (details,),
        Path("signaling/pathways") / pathway,
        "05_top_lr_pairs.png",
    )
