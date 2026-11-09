"""Optional analyses of the saved fit."""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import nnls
from scknife import run_full_pipeline
from .core import output_dir
from ._logging import _run_log


@_run_log("projection.log")
def rank_projection(fit, graph, ranks=(8, 16, 24)):
    """Fit training-only models and return reconstruction errors on the held-out cells."""
    OUT = output_dir()
    adata = fit
    G = graph
    rng = np.random.default_rng(42)
    order = rng.permutation(adata.n_obs)
    train = adata[order[:int(0.8*adata.n_obs)]].copy()
    test = adata[order[int(0.8*adata.n_obs):]].copy()
    # Remove fitted results before fitting each training-only model.
    train.uns.clear()
    for col in [c for c in train.obs if c.startswith("H_graph_")]:
        del train.obs[col]
    selected_params = json.loads((OUT / "tuning/selected_params.json").read_text())
    penalties = {
        key: selected_params[key] for key in ["lambda_reg", "beta_W", "beta_H"]
    }
    rows = []
    projections = []
    for rank in ranks:
        model = run_full_pipeline(
            train.copy(),
            G,
            k=rank,
            **penalties,
            n_iter=60,
            tol=0,
            random_state=42,
            verbose=False,
            show_loss=False,
        )
        basis = np.asarray(model.uns["Rg"] @ model.uns["W_graph"])
        genes = model.var_names
        expression = test[:, genes].X
        expression = (
            expression.toarray()
            if hasattr(expression, "toarray")
            else np.asarray(expression)
        )
        predicted = np.column_stack([nnls(basis, x)[0] for x in expression])
        reconstructed = (basis @ predicted).T
        error = np.linalg.norm(expression - reconstructed) / np.linalg.norm(expression)
        rows.append(
            {
                "rank": rank,
                "training_cells": train.n_obs,
                "heldout_cells": test.n_obs,
                "relative_reconstruction_error": float(error),
            }
        )
        if rank == selected_params["k"]:
            projections = pd.DataFrame(
                predicted.T,
                index=test.obs_names,
                columns=[f"H_projected_{i}" for i in range(rank)],
            )
        plt.close("all")
    pd.DataFrame(rows).to_csv(OUT / "rank_sensitivity.csv", index=False)
    if len(projections):
        projections.to_csv(OUT / "heldout_projected_factors.csv")
    return pd.DataFrame(rows)


def plot_projection(errors):
    """Plot relative held-out reconstruction error by rank."""
    OUT = output_dir()
    rows = errors.to_dict("records")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(
        [r["rank"] for r in rows],
        [r["relative_reconstruction_error"] for r in rows],
        marker="o",
    )
    ax.set_xticks(errors["rank"])
    ax.set_xlabel("Number of factors")
    ax.set_ylabel("Held-out relative reconstruction error")
    ax.set_title(f"Exploratory rank sensitivity · {int(errors.heldout_cells.iloc[0]):,} held-out cells")
    fig.tight_layout()
    fig.savefig(OUT / "rank_sensitivity.png", dpi=160)
    return ax


def state_programs(fit, activity, graph, n=8):
    """Return the highest pooled-activity state programs summarized by cluster."""
    OUT = output_dir()
    G = graph

    def nodes_of_kind(kind):
        return [n for n, attrs in G.nodes(data=True) if attrs.get("kind") == kind]

    nodes = nodes_of_kind("cell_state")
    pooled = activity.loc[nodes].mean(axis=1)
    selected = pooled.sort_values(ascending=False, kind="stable").head(n).index
    means = (
        activity.loc[selected]
        .T.groupby(fit.obs["leiden_clusters"], observed=True)
        .mean()
        .T
    )
    means.to_csv(OUT / "state_activity_by_cluster.csv")
    return means


def plot_states(means):
    """Show mean state-program activity across fitted clusters."""
    OUT = output_dir()
    fig, ax = plt.subplots(figsize=(11, 6))
    im = ax.imshow(means, aspect="auto", cmap="viridis")
    ax.set_xticks(range(means.shape[1]))
    ax.set_xticklabels(means.columns)
    ax.set_yticks(range(len(means)))
    ax.set_yticklabels(means.index)
    ax.set_xlabel("H-space cluster")
    ax.set_title("Cell-state program activity")
    fig.colorbar(im, ax=ax, label="Mean model-derived activity")
    fig.tight_layout()
    fig.savefig(OUT / "state_activity_by_cluster.png", dpi=160, bbox_inches="tight")
    return ax
