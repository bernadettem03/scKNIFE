# Cell-type retention

[Previous: Run full scKNIFE](10000-setup.md) · [Next: Metabolic analysis](10000-metabolic.md)

Inspect cell-type organization in the fitted representation using two views of the same cells.

## Biological question

Compare H-space clusters with the dataset’s existing cell-type labels. The model fits expression with the knowledge graph; source labels are used for sampling and visualization.

scKNIFE does not assign cell-type labels in this tutorial.

## Inputs and code

Complete [model fitting](10000-setup.md). Load its H-space clusters, UMAP and original source cell-type labels. Run from the repository root:

## Compare the learned organization with source labels

```python
import os
from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 12, "axes.titlesize": 14,
                     "axes.labelsize": 12, "xtick.labelsize": 11,
                     "ytick.labelsize": 11, "legend.fontsize": 11,
                     "axes.spines.top": False, "axes.spines.right": False})
base = Path(os.environ.get("SCKNIFE_TUTORIAL_OUTPUT", "tutorial/outputs_10000"))
adata = sc.read_h5ad(base / "fitted_subset.h5ad")
assert adata.n_obs == 10000 and adata.obs_names.is_unique
assert adata.obs[["leiden_clusters", "cellType"]].notna().all().all()
clusters = adata.obs["leiden_clusters"].astype(str)
adata.obs["leiden_clusters"] = pd.Categorical(clusters)
out = base / "cellular_organization"
out.mkdir(exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
for ax, key, title in zip(axes, ["leiden_clusters", "cellType"],
                        ["H-space clusters", "Source cell-type annotations"]):
    sc.pl.umap(adata, color=key, ax=ax, show=False, frameon=False, title=title)
fig.tight_layout()
fig.savefig(out / "organization_umap.png", dpi=160, bbox_inches="tight")
plt.close(fig)
```

```{image} ../assets/tutorial-10000/organization_umap.png
:class: compact-scientific-plot
:align: center
:alt: Shared UMAP coordinates colored by H-space cluster and source cell type.
```

Shared UMAP coordinates colored by H-space cluster and source cell type.

## Interpretation

The panels show the same 10,000 cells on identical UMAP coordinates. One panel colors H-space clusters and the other colors the source cell types, allowing visual inspection of cell-type organization in the learned representation.

## Outputs and next step

Outputs are written to `cellular_organization/`. Next, explore [metabolic activity](10000-metabolic.md) using the same fitted model.
