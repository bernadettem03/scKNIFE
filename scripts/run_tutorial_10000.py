"""Fit the 10,000-cell sample, reconstruct activities and save the core outputs."""
import os, time, json, resource
from pathlib import Path
os.environ["MPLBACKEND"]="Agg"
os.environ.setdefault("SCKNIFE_TUTORIAL_OUTPUT", str(Path.cwd()/"tutorial/outputs_10000"))
started=time.perf_counter()
from pathlib import Path
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import networkx as nx
import scanpy as sc
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 12, "axes.titlesize": 14,
                     "axes.labelsize": 12, "xtick.labelsize": 11,
                     "ytick.labelsize": 11, "legend.fontsize": 11,
                     "axes.spines.top": False, "axes.spines.right": False})
from IPython.display import display, Image
from scknife import (run_full_pipeline, compute_activity_scores,
    run_h_graph_clustering_fixed)

ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents]
            if (p / "graphs/version_1/unified_graph.graphml").exists())
OUT = Path(os.environ.get("SCKNIFE_TUTORIAL_OUTPUT", ROOT / "tutorial/outputs_10000"))
OUT.mkdir(parents=True, exist_ok=True)
SEED = 42
G = nx.read_graphml(ROOT / "graphs/version_1/unified_graph.graphml")
# %% Sample the cohort with authentic metadata
raw = sc.read_h5ad(os.environ["SCKNIFE_COUNTS"])
meta = pd.read_csv(os.environ["SCKNIFE_METADATA"], index_col=0)
assert raw.obs_names.is_unique and meta.index.is_unique
assert raw.obs_names.isin(meta.index).all()
raw.obs = meta.loc[raw.obs_names].copy()
raw = raw[raw.obs["BC_type"].eq("TNBC")].copy()
assert raw.n_obs >= 10000
strata = ["patient_id", "timepoint", "cellType"]
sizes = raw.obs.groupby(strata, observed=True, sort=True).size()
# Reserve one cell per stratum, then distribute the remaining budget proportionally.
assert len(sizes) < 10000
remaining = sizes - 1
quotas = (10000 - len(sizes)) * remaining / remaining.sum()
allocation = 1 + np.floor(quotas).astype(int)
extra = 10000 - int(allocation.sum())
allocation.loc[(quotas - np.floor(quotas)).sort_values(ascending=False, kind="stable").head(extra).index] += 1
rng = np.random.default_rng(SEED)
selected = []
for key, group in raw.obs.groupby(strata, observed=True, sort=True):
    selected.extend(rng.choice(group.index.to_numpy(), int(allocation.loc[key]), replace=False))
audit = pd.DataFrame({"cohort_cells": sizes, "tutorial_cells": allocation})
audit["cohort_fraction"] = sizes / sizes.sum()
audit["tutorial_fraction"] = allocation / 10000
audit.to_csv(OUT / "sampling_audit.csv")
raw = raw[sorted(selected)].copy()
assert raw.n_obs == 10000 and raw.obs_names.is_unique
raw.obs.to_csv(OUT / "sample_metadata.csv")
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
pd.crosstab(raw.obs.patient_id, raw.obs.timepoint).plot.bar(ax=axes[0], stacked=True)
axes[0].set_ylabel("Selected cells"); axes[0].set_title("Donor and treatment representation")
raw.obs.cellType.value_counts().plot.barh(ax=axes[1]); axes[1].set_title("Reference cell types")
fig.tight_layout(); fig.savefig(OUT / "sample_representation.png", dpi=160); plt.close(fig)
print("Selected cells × input genes:", raw.shape, flush=True)



# %% Prepare expression

adata = raw.copy()
values = adata.X.data if sp.issparse(adata.X) else np.asarray(adata.X)
assert np.isfinite(values).all() and (values >= 0).all()
assert np.allclose(values, np.round(values)), "Supply unnormalized counts"
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
shared = adata.var_names.isin(G.nodes)
adata = adata[:, shared].copy()
means = np.asarray(adata.X.mean(axis=0)).ravel()
top = np.argsort(means, kind="stable")[-2000:]
adata = adata[:, np.sort(top)].copy()
assert adata.n_vars == 2000
print("Prepared cells × genes:", adata.shape)

# %% Fit the graph model
fit = run_full_pipeline(adata.copy(), G, k=16, lambda_reg=1,
    beta_W=0.01, beta_H=0.01, n_iter=100, tol=0, random_state=SEED, verbose=False)

plt.close("all")
plt.figure(figsize=(7, 4))
plt.plot(np.arange(1, len(fit.uns["loss_history"]) + 1), fit.uns["loss_history"])
plt.xlabel("Iteration")
plt.ylabel("Normalized optimization objective")
plt.tight_layout()
plt.savefig(OUT / "convergence.png", dpi=150)
plt.show()



# %% Cluster the learned cell representation
fit, clustering = run_h_graph_clustering_fixed(fit, pca_dims=6,
    n_neighbors=15, resolution=0.5, random_state=SEED)
assert np.isfinite(fit.uns["W_graph"]).all()
assert len(fit.uns["loss_history"]) == 100
clustering.to_csv(OUT / "clustering.csv", index=False)
display(clustering)

fit = compute_activity_scores(fit, batch_size=256)
activity = fit.uns["activity_scores"]
assert activity.shape == (G.number_of_nodes(), fit.n_obs)
assert activity.columns.equals(fit.obs_names)
def nodes_of_kind(kind):
    return [n for n, attrs in G.nodes(data=True) if attrs.get("kind") == kind]
# Retain the original supplied cell identities; clustering does not assign cell types.
assert fit.obs["cellType"].astype(str).equals(raw.obs.loc[fit.obs_names,"cellType"].astype(str))
# %% Export and verify results
activity.to_parquet(OUT / "activity_scores.parquet")
fit.obs.to_csv(OUT / "cell_annotations.csv")
fit.uns.pop("activity_scores", None)
fit.write_h5ad(OUT / "fitted_subset.h5ad")
reloaded = sc.read_h5ad(OUT / "fitted_subset.h5ad")
assert reloaded.shape == fit.shape
assert reloaded.uns["W_graph"].shape == (G.number_of_nodes(), 16)
print("Tutorial completed; results saved to the configured output directory.")
report={"cells":int(fit.n_obs),"genes":int(fit.n_vars),"elapsed_seconds":time.perf_counter()-started,"peak_rss_native_units":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"peak_rss_units":"KiB on Linux; bytes on macOS"}
assert fit.n_obs==10000
(OUT/"tutorial_run_report.json").write_text(json.dumps(report,indent=2))
print(report)

import hashlib, platform, importlib.metadata
report.update({"seed": SEED, "rank": 16, "iterations": 100, "python": platform.python_version(),
 "versions": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scanpy", "scipy", "scikit-learn"]},
 "graph_sha256": hashlib.sha256((ROOT / "graphs/version_1/unified_graph.graphml").read_bytes()).hexdigest()})
(OUT / "tutorial_run_report.json").write_text(json.dumps(report, indent=2))


