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
# %% Tune on a representative subset before the full fit
from scipy.optimize import nnls
import itertools
TUNE = OUT / 'tuning'; TUNE.mkdir(exist_ok=True)
def stratified_sample(obs, number, seed, minimum=1):
    sizes = obs.groupby(strata, observed=True, sort=True).size()
    reserved = sizes.clip(upper=minimum)
    remaining = sizes - reserved
    assert reserved.sum() <= number <= sizes.sum()
    quota = (number - reserved.sum()) * remaining / remaining.sum()
    allocation = reserved + np.floor(quota).astype(int)
    extra = number - int(allocation.sum())
    allocation.loc[(quota - np.floor(quota)).sort_values(ascending=False, kind='stable').head(extra).index] += 1
    rng = np.random.default_rng(seed); ids = []
    for key, group in obs.groupby(strata, observed=True, sort=True):
        ids.extend(rng.choice(group.index.to_numpy(), int(allocation.loc[key]), replace=False))
    return sorted(ids)

tuning_ids = stratified_sample(adata.obs, 2000, SEED, minimum=2)
tuning_obs = adata.obs.loc[tuning_ids]
# Reserve one fitting cell in each stratum; select validation cells from the remainder.
reserved_train = stratified_sample(tuning_obs, len(tuning_obs.groupby(strata, observed=True)), SEED)
validation_pool = tuning_obs.drop(index=reserved_train)
validation_ids = stratified_sample(validation_pool, 500, SEED + 1)
training_ids = sorted(set(tuning_ids) - set(validation_ids))
assert len(training_ids) == 1500 and not set(training_ids) & set(validation_ids)
split = pd.Series('not used for tuning', index=adata.obs_names, name='tuning_split')
split.loc[training_ids] = 'fit'; split.loc[validation_ids] = 'validation'
adata.obs[strata].join(split).to_csv(TUNE / 'subset_membership.csv')
counts = adata.obs[strata].join(split).groupby(strata + ['tuning_split'], observed=True).size().unstack(fill_value=0)
counts.to_csv(TUNE / 'stratum_counts.csv')
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, key in zip(axes, ['patient_id', 'cellType']):
    comparison = pd.DataFrame({'All 10,000': adata.obs[key].value_counts(normalize=True),
        'Tuning 2,000': tuning_obs[key].value_counts(normalize=True),
        'Fit 1,500': adata.obs.loc[training_ids, key].value_counts(normalize=True),
        'Validation 500': adata.obs.loc[validation_ids, key].value_counts(normalize=True)}).fillna(0)
    comparison.plot.barh(ax=ax); ax.set_xlabel('Cell fraction'); ax.set_ylabel(''); ax.set_title(key)
fig.tight_layout(); fig.savefig(TUNE / 'subset_representation.png', dpi=160, bbox_inches='tight'); plt.close(fig)

# Gene selection sees fitting cells only; validation cells do not influence it.
means = np.asarray(adata[training_ids].X.mean(axis=0)).ravel()
top = np.argsort(means, kind='stable')[-2000:]
adata = adata[:, np.sort(top)].copy()
assert adata.n_vars == 2000
pd.Series(adata.var_names, name='gene').to_csv(TUNE / 'selected_genes.csv', index=False)
tune_train = adata[training_ids].copy(); tune_validation = adata[validation_ids].copy()
validation_expression = tune_validation.X.toarray() if sp.issparse(tune_validation.X) else np.asarray(tune_validation.X)
SEEDS = (42, 43)
records = []; cache = {}
tuning_started = time.perf_counter()
def evaluate(params):
    key = tuple(params[x] for x in ['k', 'lambda_reg', 'beta_W', 'beta_H'])
    if key in cache:
        return cache[key]
    errors = []
    for seed in SEEDS:
        trial_started = time.perf_counter()
        model = run_full_pipeline(tune_train.copy(), G, **params, n_iter=100,
            tol=0, random_state=seed, verbose=False)
        assert model.var_names.equals(tune_validation.var_names)
        basis = np.asarray(model.uns['Rg'] @ model.uns['W_graph'])
        projected = np.column_stack([nnls(basis, x)[0] for x in validation_expression])
        error = float(np.linalg.norm(validation_expression - (basis @ projected).T) / np.linalg.norm(validation_expression))
        assert np.isfinite(error)
        records.append(dict(params, seed=seed, validation_error=error, seconds=time.perf_counter()-trial_started))
        pd.DataFrame(records).to_csv(TUNE / 'candidate_runs.csv', index=False)
        errors.append(error)
        print('TUNING', params, 'seed', seed, 'validation error', round(error, 6), flush=True)
        plt.close('all')
    result = dict(params, mean_error=float(np.mean(errors)), seed_sd=float(np.std(errors, ddof=1)))
    cache[key] = result
    return result

# Predeclared sequential search: rank/graph penalty, then W penalty, then H penalty.
joint = pd.DataFrame([evaluate({'k': k, 'lambda_reg': lam, 'beta_W': .01, 'beta_H': .01})
                      for k, lam in itertools.product([8, 16, 24], [.1, 1., 10.])])
# Prefer the smallest rank whose mean error is within 1% of the best tested error.
threshold = 1.01 * joint.mean_error.min()
chosen = joint[joint.mean_error <= threshold].sort_values(['k', 'mean_error', 'lambda_reg']).iloc[0]
selected_params = {key: int(chosen[key]) if key == 'k' else float(chosen[key]) for key in ['k','lambda_reg','beta_W','beta_H']}
penalty_sweeps = {}
for penalty in ['beta_W', 'beta_H']:
    table = pd.DataFrame([evaluate(dict(selected_params, **{penalty: value})) for value in [.001, .01, .1]])
    best = table.sort_values(['mean_error', penalty]).iloc[0]
    selected_params[penalty] = float(best[penalty]); penalty_sweeps[penalty] = table
    table.to_csv(TUNE / (penalty + '_sweep.csv'), index=False)
joint.to_csv(TUNE / 'rank_graph_sweep.csv', index=False)
pd.DataFrame(cache.values()).to_csv(TUNE / 'candidate_summary.csv', index=False)
(TUNE / 'selected_params.json').write_text(json.dumps(selected_params, indent=2))
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
for lam, frame in joint.groupby('lambda_reg'):
    axes[0].errorbar(frame.k, frame.mean_error, yerr=frame.seed_sd, marker='o', capsize=3, label=f'lambda={lam:g}')
axes[0].axhline(threshold, color='grey', ls=':', label='1% selection band')
axes[0].axvline(selected_params['k'], color='#9b456b', ls='--')
axes[0].set(xlabel='Factors', ylabel='Relative validation reconstruction error', title='Rank and graph penalty', xticks=[8,16,24]); axes[0].legend(fontsize=9)
for ax, name in zip(axes[1:], ['beta_W','beta_H']):
    frame=penalty_sweeps[name]
    ax.errorbar(frame[name], frame.mean_error, yerr=frame.seed_sd, marker='o', capsize=3)
    ax.axvline(selected_params[name], color='#9b456b', ls='--'); ax.set_xscale('log')
    ax.set(xlabel=name, ylabel='Relative validation reconstruction error', title=name + ' (sequential sweep)')
    ax.set_ylim(axes[0].get_ylim())
fig.tight_layout(); fig.savefig(TUNE / 'parameter_selection.png', dpi=160, bbox_inches='tight'); plt.close(fig)
summary = {'tuning_cells':2000, 'fitting_cells':1500, 'validation_cells':500,
    'strata_in_full_sample':int(adata.obs.groupby(strata, observed=True).ngroups),
    'strata_in_tuning':int(tuning_obs.groupby(strata, observed=True).ngroups),
    'genes_selected_on_fitting_cells_only':True, 'seeds':list(SEEDS), 'iterations':100,
    'candidate_fits':len(records), 'selected_params':selected_params,
    'selected_validation_error':evaluate(selected_params)['mean_error'],
    'seconds':time.perf_counter()-tuning_started,
    'selection_rule':'Smallest rank within 1% of best joint-grid validation error; sequential penalties minimize mean error.'}
(TUNE / 'tuning_report.json').write_text(json.dumps(summary, indent=2))
print('LOCKED PARAMETERS', selected_params, flush=True)


# %% Fit the graph model
fit = run_full_pipeline(adata.copy(), G, **selected_params, n_iter=100, tol=0, random_state=SEED, verbose=False)

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
assert reloaded.uns["W_graph"].shape == (G.number_of_nodes(), selected_params["k"])
print("Tutorial completed; results saved to the configured output directory.")
report={"cells":int(fit.n_obs),"genes":int(fit.n_vars),"elapsed_seconds":time.perf_counter()-started,"peak_rss_native_units":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"peak_rss_units":"KiB on Linux; bytes on macOS"}
assert fit.n_obs==10000
(OUT/"tutorial_run_report.json").write_text(json.dumps(report,indent=2))
print(report)

import hashlib, platform, importlib.metadata
report.update({"seed": SEED, "rank": selected_params["k"], "selected_params": selected_params, "iterations": 100, "python": platform.python_version(),
 "versions": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scanpy", "scipy", "scikit-learn"]},
 "graph_sha256": hashlib.sha256((ROOT / "graphs/version_1/unified_graph.graphml").read_bytes()).hexdigest()})
(OUT / "tutorial_run_report.json").write_text(json.dumps(report, indent=2))


