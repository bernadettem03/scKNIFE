---
orphan: true
---

# Troubleshooting

## Python version is rejected

Use Python 3.9.21 in a fresh Conda environment. The initial release pins its scientific stack for reproducibility.

## Leiden is unavailable

Install `python-igraph` and `leidenalg`, restart Python, and rerun H-space clustering.

## No genes overlap the graph

Check identifiers, capitalization, species, and aliases. Expression genes and graph gene nodes must use compatible identifiers.

## No marker or state nodes are found

Inspect `activity.index` and confirm that marker/state names exactly match graph-node identifiers.

## Activity is transposed or cells do not align

Activity is expected as graph nodes by cells. Cell identifiers must match `adata.obs_names`.

## Full analysis is too slow or uses too much memory

Tune on a representative subset, run full factorization once, and persist W, H, activity parquet files, and the fitted H5AD. Reuse them downstream.

## Read the Docs does not rerun the tutorial

The site displays saved notebook outputs. To regenerate results with your inputs, run `python scripts/run_tutorial.py` after completing {doc}`tutorial-inputs`.
