# Full scKNIFE pipeline

Run the main model on the 10,000-cell sample prepared from the source cohort:

```bash
python scripts/run_tutorial_10000.py
```

The script uses the unchanged version-1 graph, 2,000 normalized/log-transformed shared genes, parameters selected on 1,500 fitting and 500 validation cells, 100 updates and seed 42. Follow {doc}`tutorial-inputs` and {doc}`quickstart` before running. The complete executable code is in {doc}`tutorial-10000`.

The pipeline maps expression genes to graph nodes, initializes NMF on overlapping genes, propagates node embeddings, optimizes graph-regularized W and H, records losses, and compares reconstruction with standard NMF.

## Save reusable outputs

For expensive runs, save W, H, node names, activity scores, and the fitted AnnData. Downstream analyses should load these outputs rather than rerun factorization. Keep cell identifiers unchanged so activity columns and `adata.obs_names` remain alignable.

## Reproducibility and objective

`random_state=42` controls graph-NMF H initialization. The initial standard NMF uses its existing fixed seed of 42. Repeat comparisons with explicit seeds and fixed software versions; numerical agreement across different hardware is not guaranteed. `tol=0` disables early stopping.

The monitored objective is `(0.5*||X-RgWH||² + 0.5*lambda*trace(W.T L W) + beta_W*sum(W) + beta_H*sum(H)) / ||X||²`, matching the multiplicative updates. Earlier versions reported different weights in the loss calculation. The corrected stopping criterion can therefore produce different iteration counts; see the changelog.

## Tutorial and advanced examples

Use the {doc}`tutorial-10000` for the main 10,000-cell workflow. All worked plots are generated within this 10,000-cell workflow.
