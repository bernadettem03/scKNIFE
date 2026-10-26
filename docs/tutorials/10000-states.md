---
orphan: true
---

# Other activity-based analyses: cell-state programs

[Tutorial overview](../tutorial-10000.md#analysis-tasks)

Compare graph-linked cell-state activity across clusters using the saved [fit and activities](10000-setup.md). This example selects the eight state nodes with the highest pooled mean activity.

## Summarize cell-state programs across H-space clusters

```python
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 12, "axes.titlesize": 14,
                     "axes.labelsize": 12, "xtick.labelsize": 11,
                     "ytick.labelsize": 11, "legend.fontsize": 11,
                     "axes.spines.top": False, "axes.spines.right": False})
ROOT = Path.cwd()
OUT = Path(os.environ.get("SCKNIFE_TUTORIAL_OUTPUT", ROOT / "tutorial/outputs_10000"))
fit = sc.read_h5ad(OUT / "fitted_subset.h5ad")
activity = pd.read_parquet(OUT / "activity_scores.parquet")
assert fit.n_obs == 10000 and activity.columns.equals(fit.obs_names)
assert fit.obs["cellType"].notna().all()
G = nx.read_graphml(ROOT / "graphs/version_1/unified_graph.graphml")
def nodes_of_kind(kind):
    return [n for n, attrs in G.nodes(data=True) if attrs.get("kind") == kind]

nodes = nodes_of_kind("cell_state")
pooled = activity.loc[nodes].mean(axis=1)
selected = pooled.sort_values(ascending=False, kind="stable").head(8).index
means = activity.loc[selected].T.groupby(fit.obs["leiden_clusters"], observed=True).mean().T
means.to_csv(OUT / "state_activity_by_cluster.csv")
fig, ax = plt.subplots(figsize=(11, 6))
im = ax.imshow(means, aspect="auto", cmap="viridis")
ax.set_xticks(range(means.shape[1])); ax.set_xticklabels(means.columns)
ax.set_yticks(range(len(means))); ax.set_yticklabels(means.index)
ax.set_xlabel("H-space cluster"); ax.set_title("Cell-state program activity")
fig.colorbar(im, ax=ax, label="Mean model-derived activity")
fig.tight_layout(); fig.savefig(OUT / "state_activity_by_cluster.png", dpi=160, bbox_inches="tight"); plt.close(fig)

```

## Plot cell-state activity


```{image} ../assets/tutorial-10000/state_activity_by_cluster.png
:class: compact-scientific-plot
:align: center
:alt: Mean activity of selected state programs across H-space clusters.
```

Rows are selected state-program nodes and columns are H-space clusters. Colors show mean reconstructed activity. Compare each row across clusters to see where the program has higher activity; several programs can be active in the same population.



## Outputs

`state_activity_by_cluster.csv` and `state_activity_by_cluster.png` are saved in the configured output directory. Return to the [tutorial overview](../tutorial-10000.md) to choose another analysis.
