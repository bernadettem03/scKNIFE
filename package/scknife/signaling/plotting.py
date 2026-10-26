#!/usr/bin/env python3
"""Publication figures for Scknife activity-based LR communication."""

from __future__ import annotations

import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

from .common import (
    DATA_DIR,
    FIG_DIR,
    H_COLORS,
    PATHWAY_ORDER,
    ensure_directories,
    write_json,
)

# ─── Figure Parameters (user-editable) ──────────────────────────────────────
FIGURE_WIDTH_MM = 180
FIGURE_HEIGHT_MM = 140
WIDE_FIGURE_WIDTH_MM = 260
WIDE_FIGURE_HEIGHT_MM = 145
TEXT_SIZE_PT = 8
AXIS_TITLE_SIZE_PT = 9
PANEL_LABEL_SIZE_PT = 10
TICK_LABEL_SIZE_PT = 7
LEGEND_SIZE_PT = 7
NODE_LABEL_SIZE_PT = 7
OUTPUT_DPI = 300
TOP_LR_PAIRS = 25
TOP_NETWORK_EDGES = 40
QUANTILE_LOW = 2
QUANTILE_HIGH = 98
COLORMAP_SEQUENTIAL = "viridis"
COLORMAP_DIVERGING = "RdBu_r"
EDGE_ALPHA = 0.55
NODE_SIZE = 900
# ─────────────────────────────────────────────────────────────────────────────

MM_TO_INCH = 1.0 / 25.4
FIGSIZE = (FIGURE_WIDTH_MM * MM_TO_INCH, FIGURE_HEIGHT_MM * MM_TO_INCH)
WIDE_FIGSIZE = (
    WIDE_FIGURE_WIDTH_MM * MM_TO_INCH,
    WIDE_FIGURE_HEIGHT_MM * MM_TO_INCH,
)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": TEXT_SIZE_PT,
    "axes.titlesize": PANEL_LABEL_SIZE_PT,
    "axes.labelsize": AXIS_TITLE_SIZE_PT,
    "xtick.labelsize": TICK_LABEL_SIZE_PT,
    "ytick.labelsize": TICK_LABEL_SIZE_PT,
    "legend.fontsize": LEGEND_SIZE_PT,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save_figure(fig: plt.Figure, output_base: Path) -> None:
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_base.with_suffix(".png"),
        dpi=OUTPUT_DPI,
        bbox_inches="tight",
        facecolor="white",
    )
    fig.savefig(
        output_base.with_suffix(".pdf"),
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def clean_axis(ax: plt.Axes) -> None:
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def finite_quantiles(values: np.ndarray, diverging: bool = False):
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (-1.0, 1.0) if diverging else (0.0, 1.0)
    if diverging:
        vmax = float(np.percentile(np.abs(finite), QUANTILE_HIGH))
        vmax = max(vmax, np.finfo(float).eps)
        return -vmax, vmax
    vmin, vmax = np.percentile(finite, [QUANTILE_LOW, QUANTILE_HIGH])
    if math.isclose(float(vmin), float(vmax)):
        vmax = float(vmin) + np.finfo(float).eps
    return float(vmin), float(vmax)


def group_colors(groups: list[str]) -> list[str]:
    palette = plt.get_cmap("tab20")
    fallback = {group: palette(index % 20) for index, group in enumerate(groups)}
    return [H_COLORS.get(group, fallback[group]) for group in groups]


def plot_sender_receiver(comm: pd.DataFrame, outdir: Path, prefix: str = "01") -> None:
    outgoing = comm.sum(axis=1)
    incoming = comm.sum(axis=0)
    colors = group_colors(list(comm.index))
    fig, axes = plt.subplots(1, 2, figsize=WIDE_FIGSIZE)
    axes[0].barh(comm.index, outgoing, color=colors)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Total outgoing weighted score")
    axes[0].set_title("A  Sender strength", loc="left", fontweight="bold")
    axes[1].barh(comm.columns, incoming, color=colors)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Total incoming weighted score")
    axes[1].set_title("B  Receiver strength", loc="left", fontweight="bold")
    for ax in axes:
        clean_axis(ax)
    fig.tight_layout()
    save_figure(fig, outdir / f"{prefix}_sender_receiver")


def plot_heatmaps(
    comm: pd.DataFrame,
    fdr: pd.DataFrame | None,
    outdir: Path,
    prefix: str = "03",
    title: str = "Activity-based communication",
) -> None:
    if fdr is None:
        fig, ax = plt.subplots(figsize=FIGSIZE)
        axes = [ax]
    else:
        fig, axes = plt.subplots(1, 2, figsize=WIDE_FIGSIZE)
    vmin, vmax = finite_quantiles(comm.to_numpy())
    sns.heatmap(
        comm,
        cmap=COLORMAP_SEQUENTIAL,
        vmin=vmin,
        vmax=vmax,
        square=True,
        linewidths=0.25,
        cbar_kws={"label": "Weighted score"},
        ax=axes[0],
    )
    axes[0].set_xlabel("Receiver")
    axes[0].set_ylabel("Sender")
    axes[0].set_title(f"A  {title}", loc="left", fontweight="bold")
    if fdr is not None:
        sns.heatmap(
            -np.log10(fdr.clip(lower=1e-12)),
            cmap=COLORMAP_SEQUENTIAL,
            square=True,
            linewidths=0.25,
            cbar_kws={"label": "−log10(FDR)"},
            ax=axes[1],
        )
        axes[1].set_xlabel("Receiver")
        axes[1].set_ylabel("Sender")
        axes[1].set_title("B  Permutation significance", loc="left", fontweight="bold")
    fig.tight_layout()
    save_figure(fig, outdir / f"{prefix}_communication_heatmap")


def plot_network(
    comm: pd.DataFrame,
    outdir: Path,
    prefix: str = "06",
    title: str = "H-cluster communication network",
) -> None:
    graph = nx.DiGraph()
    stacked = comm.stack().sort_values(ascending=False)
    stacked = stacked.loc[stacked > 0].head(TOP_NETWORK_EDGES)
    for (source, target), weight in stacked.items():
        graph.add_edge(source, target, weight=float(weight))
    groups = list(comm.index)
    graph.add_nodes_from(groups)
    pos = nx.circular_layout(graph)
    weights = np.array([data["weight"] for _, _, data in graph.edges(data=True)])
    if weights.size:
        lo, hi = finite_quantiles(weights)
        denom = max(hi - lo, np.finfo(float).eps)
        widths = 0.5 + 4.5 * np.clip((weights - lo) / denom, 0, 1)
    else:
        widths = []
    fig, ax = plt.subplots(figsize=FIGSIZE)
    nx.draw_networkx_nodes(
        graph,
        pos,
        node_color=group_colors(groups),
        node_size=NODE_SIZE,
        edgecolors="white",
        linewidths=1.0,
        ax=ax,
    )
    nx.draw_networkx_labels(
        graph, pos, font_size=NODE_LABEL_SIZE_PT, font_weight="bold", ax=ax
    )
    nx.draw_networkx_edges(
        graph,
        pos,
        width=widths,
        alpha=EDGE_ALPHA,
        arrows=True,
        arrowsize=10,
        connectionstyle="arc3,rad=0.08",
        edge_color="#555555",
        ax=ax,
    )
    ax.set_title(f"A  {title}", loc="left", fontweight="bold")
    ax.axis("off")
    fig.tight_layout()
    save_figure(fig, outdir / f"{prefix}_hierarchical_network")


def plot_flow(comm: pd.DataFrame, outdir: Path, prefix: str = "02") -> None:
    groups = list(comm.index)
    outgoing = comm.sum(axis=1)
    incoming = comm.sum(axis=0)
    source_y = np.linspace(0.05, 0.95, len(groups))
    target_y = np.linspace(0.05, 0.95, len(groups))
    top = comm.stack().sort_values(ascending=False)
    top = top.loc[top > 0].head(TOP_NETWORK_EDGES)
    max_weight = max(float(top.max()) if len(top) else 0.0, np.finfo(float).eps)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    for idx, group in enumerate(groups):
        color = H_COLORS.get(group, "#777777")
        ax.scatter(0.1, source_y[idx], s=NODE_SIZE / 3, color=color, zorder=3)
        ax.scatter(0.9, target_y[idx], s=NODE_SIZE / 3, color=color, zorder=3)
        ax.text(0.06, source_y[idx], group, ha="right", va="center")
        ax.text(0.94, target_y[idx], group, ha="left", va="center")
    lookup = {group: idx for idx, group in enumerate(groups)}
    for (source, target), weight in top.items():
        arrow = FancyArrowPatch(
            (0.12, source_y[lookup[source]]),
            (0.88, target_y[lookup[target]]),
            arrowstyle="-|>",
            connectionstyle="arc3,rad=0.08",
            linewidth=0.3 + 3.0 * float(weight) / max_weight,
            color=H_COLORS.get(source, "#777777"),
            alpha=EDGE_ALPHA,
        )
        ax.add_patch(arrow)
    ax.set_title("A  Dominant signal flow", loc="left", fontweight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.tight_layout()
    save_figure(fig, outdir / f"{prefix}_sankey_flow")


def plot_bubble(
    details: pd.DataFrame,
    groups: list[str],
    outdir: Path,
    prefix: str = "04",
) -> None:
    pair_summary = details.groupby(["source", "target"], as_index=False).agg(
        total_score=("weighted_score", "sum"),
        n_interactions=("weighted_score", "size"),
    )
    max_count = max(pair_summary["n_interactions"].max(), 1)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    xmap = {group: idx for idx, group in enumerate(groups)}
    ymap = {group: idx for idx, group in enumerate(groups)}
    vmin, vmax = finite_quantiles(pair_summary["total_score"].to_numpy())
    scatter = ax.scatter(
        pair_summary["source"].map(xmap),
        pair_summary["target"].map(ymap),
        s=20 + 180 * pair_summary["n_interactions"] / max_count,
        c=pair_summary["total_score"],
        cmap=COLORMAP_SEQUENTIAL,
        vmin=vmin,
        vmax=vmax,
        alpha=0.8,
        edgecolors="white",
        linewidths=0.3,
    )
    ax.set_xticks(range(len(groups)), groups, rotation=45, ha="right")
    ax.set_yticks(range(len(groups)), groups)
    ax.set_xlabel("Sender")
    ax.set_ylabel("Receiver")
    ax.set_title("A  LR communication bubble matrix", loc="left", fontweight="bold")
    fig.colorbar(scatter, ax=ax, label="Total weighted score")
    clean_axis(ax)
    fig.tight_layout()
    save_figure(fig, outdir / f"{prefix}_bubble_matrix")


def plot_top_lr(details: pd.DataFrame, outdir: Path, prefix: str = "05") -> None:
    top = (
        details.groupby(["ligand", "receptor"], as_index=False)["weighted_score"]
        .sum()
        .nlargest(TOP_LR_PAIRS, "weighted_score")
        .sort_values("weighted_score")
    )
    fig, ax = plt.subplots(figsize=FIGSIZE)
    labels = top["ligand"] + " → " + top["receptor"]
    colors = plt.get_cmap(COLORMAP_SEQUENTIAL)(
        np.linspace(0.2, 0.9, max(len(top), 1))
    )
    ax.barh(labels, top["weighted_score"], color=colors)
    ax.set_xlabel("Summed weighted score")
    ax.set_title("A  Top ligand–receptor pairs", loc="left", fontweight="bold")
    clean_axis(ax)
    fig.tight_layout()
    save_figure(fig, outdir / f"{prefix}_top_lr_pairs")


def plot_centrality(centrality: pd.DataFrame, outdir: Path, prefix: str = "07") -> None:
    centrality = centrality.set_index("bgnmf_cell_type_highest_as")
    fields = ["out_degree", "in_degree", "betweenness", "information"]
    fig, axes = plt.subplots(2, 2, figsize=WIDE_FIGSIZE)
    for panel, (ax, field) in enumerate(zip(axes.flat, fields)):
        ax.bar(
            centrality.index,
            centrality[field],
            color=group_colors(list(centrality.index)),
        )
        ax.tick_params(axis="x", rotation=45)
        ax.set_ylabel(field.replace("_", " ").title())
        ax.set_title(
            f"{chr(65 + panel)}  {field.replace('_', ' ').title()}",
            loc="left",
            fontweight="bold",
        )
        clean_axis(ax)
    fig.tight_layout()
    save_figure(fig, outdir / f"{prefix}_network_centrality")


def plot_condition_difference(
    diff: pd.DataFrame,
    outdir: Path,
    output_name: str = "08_responder_vs_nonresponder",
    title: str = "A  Descriptive signaling difference",
    colorbar_label: str = "Responder − Non-Responder",
) -> None:
    vmin, vmax = finite_quantiles(diff.to_numpy(), diverging=True)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    sns.heatmap(
        diff,
        cmap=COLORMAP_DIVERGING,
        vmin=vmin,
        vmax=vmax,
        center=0,
        square=True,
        linewidths=0.25,
        cbar_kws={"label": colorbar_label},
        ax=ax,
    )
    ax.set_xlabel("Receiver")
    ax.set_ylabel("Sender")
    ax.set_title(
        title,
        loc="left",
        fontweight="bold",
    )
    fig.tight_layout()
    save_figure(fig, outdir / output_name)


def main() -> None:
    ensure_directories()
    required = [
        DATA_DIR / "communication_matrix.csv",
        DATA_DIR / "communication_matrix_fdr.csv",
        DATA_DIR / "interaction_details.csv.gz",
        DATA_DIR / "centrality_analysis.csv",
        DATA_DIR / "dysregulated_signaling_responder_vs_nonresponder.csv",
        DATA_DIR / "dysregulated_signaling_on_vs_pre.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing computation outputs: {missing}")

    comm = pd.read_csv(required[0], index_col=0)
    fdr = pd.read_csv(required[1], index_col=0)
    details = pd.read_csv(required[2])
    centrality = pd.read_csv(required[3])
    difference = pd.read_csv(required[4], index_col=0)
    surgery_vs_bx = pd.read_csv(required[5], index_col=0)
    groups = list(comm.index)

    general = FIG_DIR / "general"
    plot_sender_receiver(comm, general)
    plot_flow(comm, general)
    plot_heatmaps(comm, fdr, general)
    plot_bubble(details, groups, general)
    plot_top_lr(details, general)
    plot_network(comm, general)
    plot_centrality(centrality, general)
    plot_condition_difference(difference, general)
    plot_condition_difference(
        surgery_vs_bx, general,
        output_name="09_on_vs_pre",
        title="A  Descriptive signaling difference: On-treatment vs Pre-treatment",
        colorbar_label="On − Pre",
    )

    pathway_inventory = {}
    selected_pathway = os.environ.get("SCKNIFE_SELECTED_PATHWAY", "VEGF")
    if selected_pathway not in PATHWAY_ORDER:
        raise ValueError(f"SCKNIFE_SELECTED_PATHWAY must be one of {PATHWAY_ORDER}")
    for pathway in [selected_pathway]:
        frame = details.loc[details["pathway_name"].astype(str) == pathway].copy()
        if frame.empty:
            pathway_inventory[pathway] = {"records": 0, "status": "not detected"}
            continue
        path_comm = frame.pivot_table(
            index="source",
            columns="target",
            values="weighted_score",
            aggfunc="sum",
            fill_value=0.0,
        ).reindex(index=groups, columns=groups, fill_value=0.0)
        outdir = FIG_DIR / "pathways" / pathway
        plot_sender_receiver(path_comm, outdir)
        plot_flow(path_comm, outdir)
        plot_heatmaps(
            path_comm,
            None,
            outdir,
            title=f"{pathway} communication",
        )
        plot_bubble(frame, groups, outdir)
        plot_top_lr(frame, outdir)
        plot_network(
            path_comm,
            outdir,
            title=f"{pathway} H-cluster network",
        )
        pathway_inventory[pathway] = {
            "records": int(len(frame)),
            "status": "plotted",
        }

    outputs = sorted(
        str(path.relative_to(FIG_DIR))
        for path in FIG_DIR.rglob("*")
        if path.is_file()
    )
    write_json(FIG_DIR / "figure_manifest.json", {
        "rendering_backend": "matplotlib/seaborn fallback; cnsplots unavailable",
        "output_dpi": OUTPUT_DPI,
        "groups": groups,
        "pathways": pathway_inventory,
        "files": outputs,
        "interpretation": (
            "Responder and On-vs-Pre comparisons are descriptive and "
            "cell-pooled; no patient-level statistical inference is shown."
        ),
    })
    print(f"Generated {len(outputs)} figure files")


if __name__ == "__main__":
    main()
