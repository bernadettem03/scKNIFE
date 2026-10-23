"""Run from the scKNIFE repository root. Uses 1,000 fully synthetic cells for software checks only. Biological signaling is not included."""
import os, time, json, resource
from pathlib import Path
os.environ["MPLBACKEND"]="Agg"
os.environ.setdefault("SCKNIFE_TUTORIAL_OUTPUT", str(Path.cwd()/"tutorial/outputs_focused"))
started=time.perf_counter()
from pathlib import Path
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import networkx as nx
import scanpy as sc
import matplotlib.pyplot as plt
from IPython.display import display, Image
from scknife import (run_full_pipeline, compute_activity_scores,
    run_h_graph_clustering_fixed, run_celltype_assignment,
    run_cellstate_assignment, h_projection)

ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents]
            if (p / "tutorial/data/synthetic_regression_1000.h5ad").exists())
OUT = Path(os.environ.get("SCKNIFE_TUTORIAL_OUTPUT", ROOT / "tutorial/outputs"))
OUT.mkdir(parents=True, exist_ok=True)
SEED = 42
G = nx.read_graphml(ROOT / "graphs/version_1/unified_graph.graphml")
raw = sc.read_h5ad(ROOT / "tutorial/data/synthetic_regression_1000.h5ad")
print(raw.shape, G.number_of_nodes(), G.number_of_edges())
assert raw.n_obs == 1000

adata = raw.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
shared = adata.var_names.isin(G.nodes)
adata = adata[:, shared].copy()
means = np.asarray(adata.X.mean(axis=0)).ravel()
top = np.argsort(means, kind="stable")[-1000:]
adata = adata[:, np.sort(top)].copy()
assert adata.n_vars == 1000
print("Prepared cells × genes:", adata.shape)

fit = run_full_pipeline(adata.copy(), G, k=8, lambda_reg=1,
    beta_W=0.01, beta_H=0.01, n_iter=60, tol=0, random_state=SEED, verbose=False)

plt.close("all")
plt.figure(figsize=(7, 4))
plt.plot(np.arange(1, len(fit.uns["loss_history"]) + 1), fit.uns["loss_history"])
plt.xlabel("Iteration")
plt.ylabel("Normalized optimization objective")
plt.tight_layout()
plt.savefig(OUT / "convergence.png", dpi=150)
plt.show()

fit, clustering = run_h_graph_clustering_fixed(fit, pca_dims=6,
    n_neighbors=15, resolution=0.5, random_state=SEED)
assert np.isfinite(fit.uns["W_graph"]).all()
assert len(fit.uns["loss_history"]) == 60
clustering.to_csv(OUT / "clustering.csv", index=False)
display(clustering)

fit = compute_activity_scores(fit, batch_size=256)
activity = fit.uns["activity_scores"]
assert activity.shape == (G.number_of_nodes(), fit.n_obs)
assert activity.columns.equals(fit.obs_names)
def nodes_of_kind(kind):
    return [n for n, attrs in G.nodes(data=True) if attrs.get("kind") == kind]
fit, celltype_tables, _ = run_celltype_assignment(fit, activity,
    nodes_of_kind("cell_type"), methods=("highest_as",))
fit, state_table, state_scores = run_cellstate_assignment(fit, activity,
    nodes_of_kind("cell_state"), cluster_key="leiden_clusters")
celltype_tables["highest_as"].to_csv(OUT / "celltype_mapping.csv", index=False)
state_table.to_csv(OUT / "cellstate_mapping.csv", index=False)
display(celltype_tables["highest_as"])
display(state_table)

sc.pl.umap(fit, color=["leiden_clusters_str", "bgnmf_cell_type_highest_as", "synthetic_group"],
           show=False, ncols=1, title=["H-space clusters", "Candidate cell types", "Synthetic groups"])
plt.savefig(OUT / "cluster_and_celltype_umap.png", dpi=150, bbox_inches="tight")
plt.show()

def selected_cluster_activity(kind, limit=8):
    node_ids = sorted(nodes_of_kind(kind))
    means = activity.loc[node_ids].T.groupby(
        fit.obs["leiden_clusters"], observed=True).mean().T
    chosen = means.mean(axis=1).sort_values(
        ascending=False, kind="stable").head(limit).index
    return means.loc[chosen]

def plot_activity_heatmap(values, title, filename):
    if values.empty:
        raise ValueError("No matching graph nodes available")
    import textwrap
    fig, ax = plt.subplots(figsize=(10, max(5, 0.65 * len(values))))
    im = ax.imshow(values.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(values.shape[1]))
    ax.set_xticklabels(values.columns.astype(str))
    ax.set_yticks(np.arange(values.shape[0]))
    ax.set_yticklabels([textwrap.fill(str(x), 42) for x in values.index], fontsize=9)
    ax.set_xlabel("Cell cluster")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="Mean model-derived activity")
    fig.tight_layout()
    fig.savefig(OUT / filename, dpi=150, bbox_inches="tight")
    plt.show()

state_means = selected_cluster_activity("cell_state")
state_means.to_csv(OUT / "selected_cell_states.csv")
plot_activity_heatmap(state_means, "Selected cell-state activity", "cell_states.png")

task_means = activity.loc[nodes_of_kind("task")].T.groupby(
    fit.obs["leiden_clusters"], observed=True).mean().T
task_means.to_csv(OUT / "metabolic_task_means.csv")
selected_tasks = selected_cluster_activity("task")
selected_tasks.to_csv(OUT / "selected_metabolic_tasks.csv")
plot_activity_heatmap(selected_tasks, "Selected metabolic-task activity", "metabolic_tasks.png")

activity.to_parquet(OUT / "activity_scores.parquet")
fit.obs.to_csv(OUT / "cell_annotations.csv")
fit.uns.pop("activity_scores", None)
fit.write_h5ad(OUT / "fitted_subset.h5ad")
reloaded = sc.read_h5ad(OUT / "fitted_subset.h5ad")
assert reloaded.shape == fit.shape
assert reloaded.uns["W_graph"].shape == (G.number_of_nodes(), 8)
print("Tutorial completed; results saved to the configured output directory.")
report={"cells":int(fit.n_obs),"genes":int(fit.n_vars),"elapsed_seconds":time.perf_counter()-started,"peak_rss_native_units":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"peak_rss_units":"KiB on Linux; bytes on macOS"}
assert fit.n_obs==1000
(OUT/"focused_run_report.json").write_text(json.dumps(report,indent=2))
print(report)
