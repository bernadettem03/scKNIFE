# Clustering and cell-type retention

[Run scKNIFE](10000-setup.md) · [Choose an analysis](../tutorial-10000.md#analysis-tasks)

Compare the saved H-space clusters with the dataset’s cell-type labels. scKNIFE does not assign new cell-type labels.

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

Both panels show the same 10,000 cells on identical UMAP coordinates, colored by cluster or source cell type.

The plot is saved in `cellular_organization/`. [Choose an analysis](../tutorial-10000.md#analysis-tasks) to explore the fitted activities.
