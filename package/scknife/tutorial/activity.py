"""MEAI, GO-node activities and exploratory within-cell-type comparisons."""

import hashlib
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ._plot_style import FIGURE_STYLE
from .core import output_dir
from scknife.tutorial_differential import (
    differential_activity,
    effect_overview,
    passes_cutoffs,
    volcano,
)


def _folder(kind):
    folder = (
        output_dir()
        / {"task": "metabolic_walkthrough", "GO_term": "go_differential"}[kind]
    )
    folder.mkdir(exist_ok=True)
    return folder


def node_scores(activity, graph, kind):
    """Return cells by graph nodes of one kind, sorted by node name."""
    nodes = sorted(
        n
        for n, d in graph.nodes(data=True)
        if d.get("kind") == kind and n in activity.index
    )
    if not nodes:
        raise ValueError(f"No scored nodes of type {kind}")
    return activity.loc[nodes].T


def calculate_meai(activity, graph):
    """Return per-cell MEAI using the pinned human task support tables."""
    BASE = output_dir()
    OUT = _folder("task")
    from urllib.request import urlopen

    TASK_TABLES = BASE / "task_support"
    TASK_TABLES.mkdir(exist_ok=True)
    TABLE_BASE = (
        "https://raw.githubusercontent.com/earmingol/scCellFie/"
        "c66c6be030a8ac871d998f17e592477b97781bf9/task_data/homo_sapiens/"
    )
    TABLE_HASHES = {
        "Task_by_Rxn.csv": "dc519e793a83a3fe98874e8ae6636b83d2b12bc3c25b42bd094dd408ba923209",
        "Task_by_Gene.csv": "79059bc0ceec00763e9fa568dc7e2d4fc8f0abc29582da29abd08910a8bf5e66",
    }
    for filename, expected in TABLE_HASHES.items():
        path = TASK_TABLES / filename
        data = (
            path.read_bytes()
            if path.exists()
            else urlopen(TABLE_BASE + filename, timeout=60).read()
        )
        assert hashlib.sha256(data).hexdigest() == expected, (
            f"Unexpected support table: {filename}"
        )
        if not path.exists():
            path.write_bytes(data)
    task_rxn = pd.read_csv(TASK_TABLES / "Task_by_Rxn.csv", index_col=0)
    task_gene = pd.read_csv(TASK_TABLES / "Task_by_Gene.csv", index_col=0)
    tasks = node_scores(activity, graph, "task").columns.tolist()

    EPS = 1e-10
    support_rows, combined = [], []
    for task in tasks:
        reactions = task_rxn.columns[task_rxn.loc[task].eq(1)].tolist()
        genes = task_gene.columns[task_gene.loc[task].eq(1)].tolist()
        assert set(reactions + genes).issubset(activity.index), (
            f"Missing activity nodes for {task}"
        )
        assert reactions and genes, f"Missing scored support for {task}"
        reaction_support = activity.loc[reactions].mean(axis=0).to_numpy()
        gene_support = activity.loc[genes].mean(axis=0).to_numpy()
        task_activity = activity.loc[task].to_numpy()
        combined.append(
            (task_activity + EPS)
            * np.sqrt((reaction_support + EPS) * (gene_support + EPS))
        )
        support_rows.append(
            {
                "task": task,
                "reactions": len(reactions),
                "genes": len(genes),
                "reaction_nodes": ";".join(reactions),
                "gene_nodes": ";".join(genes),
            }
        )
    combined = np.column_stack(combined)
    meai = pd.DataFrame(
        combined / combined.sum(axis=1, keepdims=True),
        index=activity.columns,
        columns=tasks,
    )
    assert np.isfinite(meai.to_numpy()).all() and (meai.to_numpy() >= 0).all()
    assert np.allclose(meai.sum(axis=1), 1)
    pd.DataFrame(support_rows).to_csv(OUT / "meai_support.csv", index=False)
    meai.to_parquet(OUT / "meai_scores.parquet")
    return meai


@plt.rc_context(FIGURE_STYLE)
def plot_task_umap(
    fit, scores, task="Glucose to lactate conversion", *, score_name="MEAI"
):
    """Color the fitted UMAP by one task or a sequence of tasks; return its plotting axes."""
    OUT = _folder("task")
    UMAP_TASKS = [task] if isinstance(task, str) else list(task)
    if not UMAP_TASKS or not set(UMAP_TASKS).issubset(scores.columns):
        raise ValueError("Choose tasks present in the score table")
    metabolic_scores = scores
    SCORE_NAME = score_name
    assert "X_umap" in fit.obsm, "Use the saved UMAP from the fitted tutorial model."
    xy = np.asarray(fit.obsm["X_umap"])
    assert xy.shape == (fit.n_obs, 2) and np.isfinite(xy).all()
    if not metabolic_scores.index.equals(fit.obs_names):
        raise ValueError("Score rows must match the fitted cell order.")
    fig, axes = plt.subplots(
        1,
        len(UMAP_TASKS),
        figsize=(3.25 * len(UMAP_TASKS), 3.0),
        squeeze=False,
        constrained_layout=True,
    )
    for ax, task in zip(axes.ravel(), UMAP_TASKS):
        scores = metabolic_scores[task].to_numpy()
        order = np.argsort(scores, kind="stable")
        points = ax.scatter(
            xy[order, 0],
            xy[order, 1],
            c=scores[order],
            s=4,
            cmap="viridis",
            vmin=0,
            vmax=scores.max(),
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(textwrap.fill(task, 28), fontsize=11)
        ax.set_aspect("equal")
        ax.set_axis_off()
        colorbar = fig.colorbar(points, ax=ax, fraction=0.035, pad=0.025, shrink=0.75)
        colorbar.set_label(SCORE_NAME, fontsize=10)
        colorbar.ax.tick_params(labelsize=9)
    pd.DataFrame(xy, index=fit.obs_names, columns=["UMAP1", "UMAP2"]).join(
        metabolic_scores[UMAP_TASKS]
    ).to_csv(OUT / "displayed_task_umaps.csv")
    fig.savefig(OUT / "meai_task_umaps.png", dpi=170, bbox_inches="tight")
    return axes[0, 0] if len(UMAP_TASKS) == 1 else axes[0]


_DEFAULT_SUBSET = object()


def _resolve_subset(subset, expansion):
    if subset is _DEFAULT_SUBSET:
        return {"expansion": expansion if expansion is not None else "E"}
    if expansion is not None:
        raise ValueError("Specify either subset or expansion.")
    return {} if subset is None else dict(subset)


def _subset_label(subset):
    if subset == {"expansion": "E"}:
        return "Expanded"
    if subset == {"expansion": "NE"}:
        return "Non-expanded"
    return ", ".join(f"{key}={value}" for key, value in subset.items()) or "All cells"


def _comparison_cells(fit, subset, contrast, min_cells, cell_type_key, condition_key, folder=None):
    """Align the requested metadata to the names used by the statistical kernel."""
    if cell_type_key == condition_key or len(contrast) != 2 or contrast[0] == contrast[1]:
        raise ValueError("Provide distinct grouping columns and two distinct conditions.")
    obs = fit.obs
    for key, value in subset.items():
        obs = obs.loc[obs[key].eq(value)]
    meta = obs[[cell_type_key, condition_key]].dropna().copy()
    meta.columns = ["cellType", "timepoint"]
    selected_meta = meta.loc[meta.timepoint.isin(contrast)]
    counts = pd.crosstab(selected_meta.cellType, selected_meta.timepoint)
    counts = counts.reindex(columns=list(contrast), fill_value=0)
    counts["eligible"] = counts[list(contrast)].min(axis=1).ge(min_cells)
    if folder is not None:
        counts.to_csv(folder / "condition_cell_counts.csv")
    eligible = counts.index[counts.eligible]
    test_meta = selected_meta.loc[selected_meta.cellType.isin(eligible)].copy()
    if test_meta.empty:
        raise ValueError("No cell type meets min_cells in both conditions for this subset.")
    return test_meta


@plt.rc_context(FIGURE_STYLE)
def plot_task_violin(
    fit,
    scores,
    task,
    *,
    cell_type,
    expansion=None,
    cell_type_key="cellType",
    condition_key="timepoint",
    subset=_DEFAULT_SUBSET,
    contrast=("Pre", "On"),
    score_name="MEAI",
    min_cells=20,
):
    """Plot two score distributions for one task and source cell type; return the axes."""
    subset = _resolve_subset(subset, expansion)
    test_meta = _comparison_cells(fit, subset, contrast, min_cells, cell_type_key, condition_key)
    if cell_type not in set(test_meta.cellType):
        raise ValueError("Too few cells for the selected comparison")
    CELL_TYPE, EXAMPLE_TASK, metabolic_scores, CONTRAST, SCORE_NAME = (
        cell_type,
        task,
        scores,
        contrast,
        score_name,
    )
    GROUP_NAME = _subset_label(subset)
    example_meta = test_meta.loc[test_meta.cellType.eq(CELL_TYPE)]
    values = [
        metabolic_scores.loc[
            example_meta.index[example_meta.timepoint.eq(c)], EXAMPLE_TASK
        ].to_numpy()
        for c in CONTRAST
    ]
    positions = [0.8, 1.35]
    fig, ax = plt.subplots(figsize=(3.05, 3.5))
    violins = ax.violinplot(
        values, positions=positions, widths=0.28, showmedians=False, showextrema=False
    )
    for body, color in zip(violins["bodies"], ["#327e9b", "#c86556"]):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_linewidth(0.8)
        body.set_alpha(0.72)
    for position, value in zip(positions, values):
        q1, median, q3 = np.quantile(value, [0.25, 0.5, 0.75])
        ax.vlines(position, q1, q3, color="#243344", linewidth=2.2)
        ax.scatter(
            position,
            median,
            s=20,
            c="white",
            edgecolors="#243344",
            linewidths=0.7,
            zorder=3,
        )
    ax.set_xticks(positions, list(CONTRAST))
    ax.set_xlim(0.45, 1.7)
    ax.set_ylabel(SCORE_NAME)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)
    ax.set_title(
        f"{textwrap.fill(EXAMPLE_TASK, 26)}\n{CELL_TYPE.replace('_', ' ')} · {GROUP_NAME}",
        fontsize=10.5, pad=13,
    )
    fig.tight_layout()
    fig.savefig(
        _folder("task") / "metabolic_task_distribution.png",
        dpi=170,
        bbox_inches="tight",
    )
    return ax


def compare_activity(
    fit,
    scores,
    *,
    kind,
    score_name,
    expansion=None,
    cell_type_key="cellType",
    condition_key="timepoint",
    subset=_DEFAULT_SUBSET,
    contrast=("Pre", "On"),
    min_cells=20,
    q_cutoff=0.05,
    min_d=0.3,
    min_log2fc=0.25,
):
    """Test scores within cell types. Positive effects favor the second condition.

    cell_type_key and condition_key name columns in fit.obs. Supply a mapping of
    metadata columns to values to select cells, or subset=None to include all.
    Omitting both subset and expansion retains the expanded-cell default for the tutorial dataset.
    """
    folder = _folder(kind)
    subset = _resolve_subset(subset, expansion)
    meta = _comparison_cells(fit, subset, contrast, min_cells, cell_type_key, condition_key, folder)
    result = differential_activity(
        scores, meta, contrast, min_cells, q_cutoff, min_d, min_log2fc
    )
    result["score"] = score_name
    stem = {"task": "metabolic", "GO_term": "go"}[kind]
    result.to_csv(folder / (stem + "_differential_activity.csv"), index=False)
    result.attrs.update(
        kind=kind,
        score_name=score_name,
        expansion=subset.get("expansion"),
        subset=subset,
        contrast=contrast,
        min_d=min_d,
        q_cutoff=q_cutoff,
        min_log2fc=min_log2fc,
    )
    return result


def _plot_args(results, suffix):
    kind = results.attrs["kind"]
    stem, title = (
        ("task", "Metabolic tasks · " + results.attrs["score_name"])
        if kind == "task"
        else ("go", "GO biological processes")
    )
    return (
        title,
        _folder(kind) / (stem + "_" + suffix + ".png"),
        _subset_label(results.attrs.get("subset", {"expansion": results.attrs.get("expansion")})),
    )


def _display_results(results, q_cutoff, min_d, min_log2fc):
    if q_cutoff is None and min_d is None and min_log2fc is None:
        return results
    if min_log2fc is None and "min_log2fc" not in results.attrs:
        raise ValueError("Supply min_log2fc for results saved without this cutoff.")
    q = results.attrs["q_cutoff"] if q_cutoff is None else q_cutoff
    d = results.attrs["min_d"] if min_d is None else min_d
    fc = results.attrs["min_log2fc"] if min_log2fc is None else min_log2fc
    display = results.assign(passes_cutoffs=passes_cutoffs(results, q, d, fc))
    display.attrs.update(q_cutoff=q, min_d=d, min_log2fc=fc)
    return display


def select_results(results, *, cell_type=None, q_cutoff=None, min_d=None, min_log2fc=None):
    """Return all passing rows, optionally for one tested cell type, ordered by adjusted p-value."""
    display = _display_results(results, q_cutoff, min_d, min_log2fc)
    if cell_type is not None:
        if cell_type not in set(display.cell_type):
            raise ValueError("Cell type was not tested")
        display = display.loc[display.cell_type.eq(cell_type)]
    return display.loc[display.passes_cutoffs].sort_values("adj_p_value", kind="stable").copy()


@plt.rc_context(FIGURE_STYLE)
def plot_effects(results, *, q_cutoff=None, min_d=None, min_log2fc=None):
    """Plot effects by cell type, optionally overriding the comparison's display cutoffs.

    Highlight rows with an adjusted p-value below q_cutoff and absolute
    Cohen's d and log2 fold change at or above min_d and min_log2fc. Omitted cutoffs come from compare_activity.
    Returns the Matplotlib axes and preserves the input results table.
    """
    results = _display_results(results, q_cutoff, min_d, min_log2fc)
    title, file, group = _plot_args(results, "effects")
    return effect_overview(
        results,
        title,
        file,
        group,
        results.attrs["min_d"],
        contrast=results.attrs["contrast"],
        return_axes=True,
    )


@plt.rc_context(FIGURE_STYLE)
def plot_volcano(results, *, cell_type, q_cutoff=None, min_d=None, min_log2fc=None):
    """Plot adjusted significance against effect size for a tested cell type.

    Unspecified display cutoffs use the values from compare_activity. All three
    cutoffs determine which points are colored and eligible for labeling.
    """
    title, file, group = _plot_args(results, "volcano")
    if cell_type not in set(results.cell_type):
        raise ValueError("Cell type was not tested")
    results = _display_results(results, q_cutoff, min_d, min_log2fc)
    return volcano(
        results,
        title,
        file,
        cell_type,
        group,
        results.attrs["q_cutoff"],
        results.attrs["min_d"],
        contrast=results.attrs["contrast"],
        return_axes=True,
    )
