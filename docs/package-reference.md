# scKNIFE package reference

## Core workflow

- `run_full_pipeline(adata, G, ...)` — fit graph-regularized NMF.
- `compute_activity_scores(adata)` — reconstruct graph-node activities.
- `run_tuning(adata, G, ...)` — tune dataset-specific parameters.

## Graphs

Use `construct_unified_graph()` with Version 1, Version 2, individual supplied GraphML files, or custom NetworkX graphs. Keep the graph fixed between tuning and final fitting.

## Downstream APIs

- `run_h_graph_clustering(adata, ...)` — select H-space Leiden clustering by silhouette and compute UMAP.
- `run_h_graph_clustering_fixed(adata, ...)` — run one specified H → PCA → neighbors → Leiden → UMAP configuration without grid search.
- `run_scanpy_KNIFE_DE(...)` — differential analysis of activity features.
- `run_activity_lr(...)` — compute activity-based ligand-receptor tables with input-fingerprinted, resumable permutation checkpoints.

These signaling scores are activity proxies, not CellChat communication probabilities. Use {doc}`tutorial-10000` for the descriptive 10,000-cell scoring example; see {doc}`signaling` for the separate metadata-based API and permutation settings.
