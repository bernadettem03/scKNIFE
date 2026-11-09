"""Inputs, tuning and fitted outputs for the tutorial dataset."""

from pathlib import Path
from ._logging import _run_log
import os
import itertools
import json
import time
import numpy as np
import pandas as pd
import networkx as nx
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.optimize import nnls
from scknife import run_full_pipeline, compute_activity_scores, construct_unified_graph


def repository_root():
    """Locate the tutorial checkout, including when Jupyter starts in tutorial/."""
    for directory in (Path.cwd(), *Path.cwd().parents):
        if (directory / "graphs/version_1/unified_graph.graphml").is_file():
            return directory
    raise FileNotFoundError("Start the tutorial from the scKNIFE checkout.")


def output_dir():
    """Return the tutorial output directory, creating it on first use."""
    configured = os.environ.get("SCKNIFE_TUTORIAL_OUTPUT")
    folder = (
        Path(configured) if configured else repository_root() / "tutorial/outputs_36721"
    )
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def load_graph(path=None):
    """Load an explicit graph, or reopen the graph used by the saved tutorial fit."""
    if path is None:
        fitted = output_dir() / "fitted_graph.graphml"
        path = fitted if fitted.is_file() else "graphs/version_1/unified_graph.graphml"
    path = Path(path)
    return nx.read_graphml(path if path.is_file() else repository_root() / path)


def build_graph(collections):
    """Combine named individual collections using the scKNIFE graph builder."""
    directory = repository_root() / "graphs/individual"
    choices = {p.parent.name: p for p in directory.glob("*/*.graphml")}
    if (
        isinstance(collections, str)
        or not collections
        or not set(collections) <= choices.keys()
    ):
        raise ValueError(f"Choose one or more collections from {sorted(choices)}")
    return construct_unified_graph(
        {name: nx.read_graphml(choices[name]) for name in collections}
    )


def marker_graph(csv_path, *, source):
    """Build a cell-type marker graph from a tissue-filtered gene,cell_type CSV."""
    data = pd.read_csv(csv_path, dtype=str, usecols=["gene", "cell_type"])
    data = (
        data.apply(lambda column: column.str.strip())
        .replace("", pd.NA)
        .dropna()
        .drop_duplicates()
    )
    if data.empty:
        raise ValueError("Provide at least one gene and cell_type pair")
    graph = nx.Graph()
    for gene, label in data[["gene", "cell_type"]].itertuples(index=False, name=None):
        cell_type = f"cell_type:{label}"
        graph.add_node(gene, kind="gene")
        graph.add_node(cell_type, kind="cell_type", label=label)
        graph.add_edge(gene, cell_type, source=source, type="gene-celltype")
    return graph


def _raw_counts(path):
    """Read the original sparse counts RDS, retaining its cell and gene identifiers."""
    from tempfile import TemporaryDirectory
    from subprocess import run, DEVNULL
    from shutil import which
    from scipy.io import mmread
    from anndata import AnnData

    if path.suffix.lower() == ".h5ad":
        return sc.read_h5ad(path)
    if path.suffix.lower() != ".rds":
        raise ValueError("Expected a raw sparse RDS or unnormalized H5AD count matrix")
    if not path.is_file():
        raise FileNotFoundError(path)
    if not which("Rscript"):
        raise RuntimeError(
            "Reading the raw RDS requires Rscript and the R Matrix package."
        )
    with TemporaryDirectory(prefix="scknife-counts-") as tmp:
        folder = Path(tmp)
        run(
            [
                "Rscript",
                "--vanilla",
                str(repository_root() / "scripts/prepare_bassez_counts.R"),
                str(path),
                tmp,
            ],
            check=True,
            stdout=DEVNULL,
        )
        matrix = mmread(folder / "counts.mtx").T.tocsr()
        cells = (folder / "cells.txt").read_text().splitlines()
        genes = (folder / "genes.txt").read_text().splitlines()
    if len(set(genes)) != len(genes) or len(set(cells)) != len(cells):
        raise ValueError("Counts must have unique gene and cell identifiers")
    return AnnData(matrix, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))


def load_fit():
    """Open the saved model, including its cluster embedding once computed."""
    return sc.read_h5ad(output_dir() / "fitted_subset.h5ad")


def load_activities(fit):
    """Read the node-by-cell scores, checking identity and column order."""
    scores = pd.read_parquet(output_dir() / "activity_scores.parquet")
    if not scores.columns.equals(fit.obs_names):
        raise ValueError("Activity columns must match the fitted cells in order.")
    return scores


def save_fit(fit):
    """Save graph activities separately from the model to avoid duplicating them."""
    folder = output_dir()
    if "activity_scores" in fit.uns:
        scores = fit.uns["activity_scores"]
        assert scores.columns.equals(fit.obs_names)
        scores.to_parquet(folder / "activity_scores.parquet")
        del fit.uns["activity_scores"]
    fit.obs.to_csv(folder / "cell_annotations.csv")
    fit.write_h5ad(folder / "fitted_subset.h5ad")


def figure(relative_path, width=680):
    """Display a generated PNG in a notebook (no figure is embedded in this helper)."""
    from IPython.display import Image

    return Image(filename=str(output_dir() / relative_path), width=width)


def sample_bassez(counts=None, metadata=None, n_cells=None, seed=42):
    """Align raw cohort-2 counts and metadata and sample TNBC strata proportionally."""
    OUT = output_dir()
    if counts is None or metadata is None:
        inputs = repository_root() / "tutorial/inputs"
        if counts is None:
            counts = inputs / "1867-counts_cells_cohort2.rds"
        if metadata is None:
            metadata = inputs / "1871-BIOKEY_metaData_cohort2_web.csv"
    raw = _raw_counts(Path(counts))
    meta = pd.read_csv(metadata, index_col=0)
    assert raw.obs_names.is_unique and meta.index.is_unique
    assert raw.obs_names.isin(meta.index).all()
    raw.obs = meta.loc[raw.obs_names].copy()
    raw = raw[raw.obs["BC_type"].eq("TNBC")].copy()
    if n_cells is None:
        n_cells = raw.n_obs
    assert raw.n_obs >= n_cells
    strata = ["patient_id", "timepoint", "cellType"]
    sizes = raw.obs.groupby(strata, observed=True, sort=True).size()
    # Reserve one cell per stratum, then distribute the remaining budget proportionally.
    assert len(sizes) < n_cells
    remaining = sizes - 1
    quotas = (n_cells - len(sizes)) * remaining / remaining.sum()
    allocation = 1 + np.floor(quotas).astype(int)
    extra = n_cells - int(allocation.sum())
    allocation.loc[
        (quotas - np.floor(quotas))
        .sort_values(ascending=False, kind="stable")
        .head(extra)
        .index
    ] += 1
    rng = np.random.default_rng(seed)
    selected = []
    for key, group in raw.obs.groupby(strata, observed=True, sort=True):
        selected.extend(
            rng.choice(group.index.to_numpy(), int(allocation.loc[key]), replace=False)
        )
    audit = pd.DataFrame({"cohort_cells": sizes, "tutorial_cells": allocation})
    audit["cohort_fraction"] = sizes / sizes.sum()
    audit["tutorial_fraction"] = allocation / n_cells
    audit.to_csv(OUT / "sampling_audit.csv")
    raw = raw[sorted(selected)].copy()
    assert raw.n_obs == n_cells and raw.obs_names.is_unique
    raw.obs.to_csv(OUT / "sample_metadata.csv")
    return raw


def plot_sample(raw):
    """Plot the sampled patient/timepoint counts and cell-type counts."""
    OUT = output_dir()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    pd.crosstab(raw.obs.patient_id, raw.obs.timepoint).plot.bar(
        ax=axes[0], stacked=True
    )
    axes[0].set_ylabel("Selected cells")
    axes[0].set_title("Donor and treatment representation")
    raw.obs.cellType.value_counts().plot.barh(ax=axes[1])
    axes[1].set_title("Reference cell types")
    fig.tight_layout()
    fig.savefig(OUT / "sample_representation.png", dpi=160)
    return axes


def prepare_tuning(raw, graph, subset=2000, validation=500, seed=42):
    """Normalize graph-shared expression and choose disjoint stratified tuning partitions."""
    strata = ["patient_id", "timepoint", "cellType"]
    OUT = output_dir()
    adata = raw.copy()
    values = adata.X.data if sp.issparse(adata.X) else np.asarray(adata.X)
    assert np.isfinite(values).all() and (values >= 0).all()
    assert np.allclose(values, np.round(values)), "Supply unnormalized counts"
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    shared = adata.var_names.isin(graph.nodes)
    adata = adata[:, shared].copy()
    # %% Tune on a representative subset before the full fit

    TUNE = OUT / "tuning"
    TUNE.mkdir(exist_ok=True)

    def stratified_sample(obs, number, seed, minimum=1):
        sizes = obs.groupby(strata, observed=True, sort=True).size()
        reserved = sizes.clip(upper=minimum)
        remaining = sizes - reserved
        assert reserved.sum() <= number <= sizes.sum()
        quota = (number - reserved.sum()) * remaining / remaining.sum()
        allocation = reserved + np.floor(quota).astype(int)
        extra = number - int(allocation.sum())
        allocation.loc[
            (quota - np.floor(quota))
            .sort_values(ascending=False, kind="stable")
            .head(extra)
            .index
        ] += 1
        rng = np.random.default_rng(seed)
        ids = []
        for key, group in obs.groupby(strata, observed=True, sort=True):
            ids.extend(
                rng.choice(
                    group.index.to_numpy(), int(allocation.loc[key]), replace=False
                )
            )
        return sorted(ids)

    tuning_ids = stratified_sample(adata.obs, subset, seed, minimum=2)
    tuning_obs = adata.obs.loc[tuning_ids]
    # Reserve one fitting cell in each stratum; select validation cells from the remainder.
    reserved_train = stratified_sample(
        tuning_obs, len(tuning_obs.groupby(strata, observed=True)), seed
    )
    validation_pool = tuning_obs.drop(index=reserved_train)
    validation_ids = stratified_sample(validation_pool, validation, seed + 1)
    training_ids = sorted(set(tuning_ids) - set(validation_ids))
    assert len(training_ids) == subset - validation and not set(training_ids) & set(
        validation_ids
    )
    split = pd.Series("not used for tuning", index=adata.obs_names, name="tuning_split")
    split.loc[training_ids] = "fit"
    split.loc[validation_ids] = "validation"
    adata.obs[strata].join(split).to_csv(TUNE / "subset_membership.csv")
    counts = (
        adata.obs[strata]
        .join(split)
        .groupby(strata + ["tuning_split"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    counts.to_csv(TUNE / "stratum_counts.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, key in zip(axes, ["patient_id", "cellType"]):
        comparison = pd.DataFrame(
            {
                f"All {len(adata):,}": adata.obs[key].value_counts(normalize=True),
                f"Tuning {subset:,}": tuning_obs[key].value_counts(normalize=True),
                f"Fit {subset - validation:,}": adata.obs.loc[
                    training_ids, key
                ].value_counts(normalize=True),
                f"Validation {validation:,}": adata.obs.loc[
                    validation_ids, key
                ].value_counts(normalize=True),
            }
        ).fillna(0)
        comparison.plot.barh(ax=ax)
        ax.set_xlabel("Cell fraction")
        ax.set_ylabel("")
        ax.set_title(key)
    fig.tight_layout()
    fig.savefig(TUNE / "subset_representation.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return adata, training_ids, validation_ids


@_run_log("tuning/fit.log")
def select_parameters(
    adata,
    graph,
    training_ids,
    validation_ids,
    *,
    ranks=(8, 16, 24),
    graph_penalties=(0.1, 1, 10),
    penalties=(0.001, 0.01, 0.1),
    n_genes=2000,
    seeds=(42, 43),
):
    """Tune by held-out reconstruction; return the gene-selected expression and chosen parameters."""
    G = graph
    TUNE = output_dir() / "tuning"
    strata = ["patient_id", "timepoint", "cellType"]
    # Gene selection sees fitting cells only; validation cells do not influence it.
    means = np.asarray(adata[training_ids].X.mean(axis=0)).ravel()
    top = np.argsort(means, kind="stable")[-n_genes:]
    adata = adata[:, np.sort(top)].copy()
    assert adata.n_vars == n_genes
    pd.Series(adata.var_names, name="gene").to_csv(
        TUNE / "selected_genes.csv", index=False
    )
    tune_train = adata[training_ids].copy()
    tune_validation = adata[validation_ids].copy()
    validation_expression = (
        tune_validation.X.toarray()
        if sp.issparse(tune_validation.X)
        else np.asarray(tune_validation.X)
    )
    SEEDS = seeds
    records = []
    cache = {}
    tuning_started = time.perf_counter()

    def evaluate(params):
        key = tuple(params[x] for x in ["k", "lambda_reg", "beta_W", "beta_H"])
        if key in cache:
            return cache[key]
        errors = []
        for seed in SEEDS:
            trial_started = time.perf_counter()
            model = run_full_pipeline(
                tune_train.copy(),
                G,
                **params,
                n_iter=100,
                tol=0,
                random_state=seed,
                verbose=False,
                show_loss=False,
            )
            assert model.var_names.equals(tune_validation.var_names)
            basis = np.asarray(model.uns["Rg"] @ model.uns["W_graph"])
            projected = np.column_stack(
                [nnls(basis, x)[0] for x in validation_expression]
            )
            error = float(
                np.linalg.norm(validation_expression - (basis @ projected).T)
                / np.linalg.norm(validation_expression)
            )
            assert np.isfinite(error)
            records.append(
                dict(
                    params,
                    seed=seed,
                    validation_error=error,
                    seconds=time.perf_counter() - trial_started,
                )
            )
            pd.DataFrame(records).to_csv(TUNE / "candidate_runs.csv", index=False)
            errors.append(error)
            print(
                "TUNING",
                params,
                "seed",
                seed,
                "validation error",
                round(error, 6),
                flush=True,
            )
            plt.close("all")
        result = dict(
            params,
            mean_error=float(np.mean(errors)),
            seed_sd=float(np.std(errors, ddof=1)),
        )
        cache[key] = result
        return result

    # Predeclared sequential search: rank/graph penalty, then W penalty, then H penalty.
    joint = pd.DataFrame(
        [
            evaluate({"k": k, "lambda_reg": lam, "beta_W": 0.01, "beta_H": 0.01})
            for k, lam in itertools.product(ranks, graph_penalties)
        ]
    )
    # Prefer the smallest rank whose mean error is within 1% of the best tested error.
    threshold = 1.01 * joint.mean_error.min()
    chosen = (
        joint[joint.mean_error <= threshold]
        .sort_values(["k", "mean_error", "lambda_reg"])
        .iloc[0]
    )
    selected_params = {
        key: int(chosen[key]) if key == "k" else float(chosen[key])
        for key in ["k", "lambda_reg", "beta_W", "beta_H"]
    }
    penalty_sweeps = {}
    for penalty in ["beta_W", "beta_H"]:
        table = pd.DataFrame(
            [evaluate(dict(selected_params, **{penalty: value})) for value in penalties]
        )
        best = table.sort_values(["mean_error", penalty]).iloc[0]
        selected_params[penalty] = float(best[penalty])
        penalty_sweeps[penalty] = table
        table.to_csv(TUNE / (penalty + "_sweep.csv"), index=False)
    joint.to_csv(TUNE / "rank_graph_sweep.csv", index=False)
    pd.DataFrame(cache.values()).to_csv(TUNE / "candidate_summary.csv", index=False)
    (TUNE / "selected_params.json").write_text(json.dumps(selected_params, indent=2))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for lam, frame in joint.groupby("lambda_reg"):
        axes[0].errorbar(
            frame.k,
            frame.mean_error,
            yerr=frame.seed_sd,
            marker="o",
            capsize=3,
            label=f"lambda={lam:g}",
        )
    axes[0].axhline(threshold, color="grey", ls=":", label="1% selection band")
    axes[0].axvline(selected_params["k"], color="#9b456b", ls="--")
    axes[0].set(
        xlabel="Factors",
        ylabel="Relative validation reconstruction error",
        title="Rank and graph penalty",
        xticks=ranks,
    )
    axes[0].legend(fontsize=9)
    for ax, name in zip(axes[1:], ["beta_W", "beta_H"]):
        frame = penalty_sweeps[name]
        ax.errorbar(
            frame[name], frame.mean_error, yerr=frame.seed_sd, marker="o", capsize=3
        )
        ax.axvline(selected_params[name], color="#9b456b", ls="--")
        ax.set_xscale("log")
        ax.set(
            xlabel=name,
            ylabel="Relative validation reconstruction error",
            title=name + " (sequential sweep)",
        )
        ax.set_ylim(axes[0].get_ylim())
    fig.tight_layout()
    fig.savefig(TUNE / "parameter_selection.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    summary = {
        "tuning_cells": len(training_ids) + len(validation_ids),
        "fitting_cells": len(training_ids),
        "validation_cells": len(validation_ids),
        "strata_in_full_sample": int(adata.obs.groupby(strata, observed=True).ngroups),
        "strata_in_tuning": int(
            adata.obs.loc[training_ids + validation_ids]
            .groupby(strata, observed=True)
            .ngroups
        ),
        "genes_selected_on_fitting_cells_only": True,
        "seeds": list(SEEDS),
        "iterations": 100,
        "candidate_fits": len(records),
        "selected_params": selected_params,
        "selected_validation_error": evaluate(selected_params)["mean_error"],
        "seconds": time.perf_counter() - tuning_started,
        "selection_rule": "Smallest rank within 1% of best joint-grid validation error; sequential penalties minimize mean error.",
    }
    (TUNE / "tuning_report.json").write_text(json.dumps(summary, indent=2))
    print("LOCKED PARAMETERS", selected_params, flush=True)
    return adata, selected_params


def plot_convergence(fit):
    """Plot the optimization objective over the full fit."""
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(np.arange(1, len(fit.uns["loss_history"]) + 1), fit.uns["loss_history"])
    ax.set(xlabel="Iteration", ylabel="Normalized optimization objective")
    fig.tight_layout()
    fig.savefig(output_dir() / "convergence.png", dpi=150)
    return ax


def plot_clusters(fit):
    """Color the existing H-space UMAP by its fitted Leiden clusters."""
    fig, ax = plt.subplots(figsize=(6, 5))
    from anndata import AnnData

    # Cluster IDs denote categories; integer coloring would imply a continuous scale.
    labels = fit.obs["leiden_clusters"].astype(str)
    plot_data = AnnData(
        obs=pd.DataFrame({"cluster": pd.Categorical(labels)}, index=fit.obs_names)
    )
    plot_data.obsm["X_umap"] = fit.obsm["X_umap"]
    sc.pl.umap(
        plot_data,
        color="cluster",
        ax=ax,
        show=False,
        frameon=False,
        title="H-space clusters",
    )
    folder = output_dir() / "cellular_organization"
    folder.mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(folder / "organization_umap.png", dpi=160, bbox_inches="tight")
    return ax


def optimize_hyperparameters(
    sample,
    graph,
    *,
    ranks=(8, 16, 24),
    graph_penalties=(0.1, 1, 10),
    penalties=(0.001, 0.01, 0.1),
    subset=2000,
    validation=500,
    seed=42,
):
    """Prepare graph-shared expression and select one set of model parameters."""
    adata, train, held_out = prepare_tuning(sample, graph, subset, validation, seed)
    return select_parameters(
        adata,
        graph,
        train,
        held_out,
        ranks=ranks,
        graph_penalties=graph_penalties,
        penalties=penalties,
        seeds=(seed, seed + 1),
    )


@_run_log("full_fit.log")
def fit_and_score(adata, graph, params, *, n_iter=100, seed=42):
    """Fit the full sample with the selected settings; save the fit, scores and chosen graph."""
    fit = run_full_pipeline(
        adata.copy(),
        graph,
        **params,
        n_iter=n_iter,
        tol=0,
        random_state=seed,
        verbose=False,
        show_loss=False,
    )
    fit = compute_activity_scores(fit, batch_size=256)
    save_fit(fit)
    nx.write_graphml(graph, output_dir() / "fitted_graph.graphml")
    return fit
