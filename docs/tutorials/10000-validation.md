# Optional: rank sensitivity and projection

[Tutorial overview](../tutorial-10000.md) · [Previous: Cell–cell signaling](10000-signaling.md)

**Input:** the saved object from [model fitting](10000-setup.md). Use the command below to generate the projection results, or use the Python example to adjust the settings.

```bash
python scripts/run_tutorial_10000_validation.py
```

Split the selected sample into **8,000 training and 2,000 held-out cells**, seed 42. Fit ranks 8, 16 and 24 on training cells and project held-out expression onto each fixed `Rg @ W` basis using nonnegative least squares. The primary biological analysis remains the separate tuned fit on all 10,000 cells.

The pre-fit tuning settings supply the W/H and graph penalties. This optional split illustrates projection within the same cohort; it is not a second independent evaluation of the tuning choices. Projected H values are exported for the selected model rank.

## Rank sensitivity and projection code

```python
"""Exploratory rank sensitivity and held-out projection within the tutorial sample."""
import os,time,json,resource
from pathlib import Path
os.environ['MPLBACKEND']='Agg'
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 12, "axes.titlesize": 14,
                     "axes.labelsize": 12, "xtick.labelsize": 11,
                     "ytick.labelsize": 11, "legend.fontsize": 11,
                     "axes.spines.top": False, "axes.spines.right": False})
from scipy.optimize import nnls
from scknife import run_full_pipeline
start=time.perf_counter()
ROOT=Path.cwd(); OUT=Path(os.environ.get('SCKNIFE_TUTORIAL_OUTPUT', ROOT/'tutorial/outputs_10000'))
adata=sc.read_h5ad(OUT/'fitted_subset.h5ad')
G=nx.read_graphml(ROOT/'graphs/version_1/unified_graph.graphml')
rng=np.random.default_rng(42); order=rng.permutation(adata.n_obs)
train=adata[order[:8000]].copy(); test=adata[order[8000:]].copy()
# Remove fitted results before fitting each training-only model.
train.uns.clear()
for col in [c for c in train.obs if c.startswith('H_graph_')]: del train.obs[col]
selected_params=json.loads((OUT/'tuning/selected_params.json').read_text())
penalties={key:selected_params[key] for key in ['lambda_reg','beta_W','beta_H']}
rows=[]; projections=[]
for rank in (8,16,24):
    model=run_full_pipeline(train.copy(),G,k=rank,**penalties,n_iter=60,tol=0,random_state=42,verbose=False)
    basis=np.asarray(model.uns['Rg'] @ model.uns['W_graph'])
    genes=model.var_names
    expression=test[:,genes].X
    expression=expression.toarray() if hasattr(expression,'toarray') else np.asarray(expression)
    predicted=np.column_stack([nnls(basis,x)[0] for x in expression])
    reconstructed=(basis@predicted).T
    error=np.linalg.norm(expression-reconstructed)/np.linalg.norm(expression)
    rows.append({'rank':rank,'training_cells':8000,'heldout_cells':2000,'relative_reconstruction_error':float(error)})
    if rank==selected_params["k"]:
        projections=pd.DataFrame(predicted.T,index=test.obs_names,columns=[f'H_projected_{i}' for i in range(rank)])
    plt.close('all')
pd.DataFrame(rows).to_csv(OUT/'rank_sensitivity.csv',index=False)
projections.to_csv(OUT/'heldout_projected_factors.csv')
fig,ax=plt.subplots(figsize=(6,4));ax.plot([r['rank'] for r in rows],[r['relative_reconstruction_error'] for r in rows],marker='o')
ax.set_xticks([8,16,24]);ax.set_xlabel('Number of factors');ax.set_ylabel('Held-out relative reconstruction error')
ax.set_title('Exploratory rank sensitivity · 2,000 held-out cells')
fig.tight_layout();fig.savefig(OUT/'rank_sensitivity.png',dpi=160);plt.close(fig)
report={'elapsed_seconds':time.perf_counter()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'training_cells':8000,'heldout_cells':2000,'gene_selection':'tuning fitting subset; overlaps this exploratory split', 'selected_params':selected_params,'purpose':'exploratory reconstruction diagnostic, not unbiased generalization or biological rank selection'}
(OUT/'validation_report.json').write_text(json.dumps(report,indent=2))
print(rows,report)
```

```{image} ../assets/tutorial-10000/rank_sensitivity.png
:class: compact-scientific-plot
:align: center
:alt: Relative held-out reconstruction error across the three tested ranks.
```

Relative held-out reconstruction error across the three tested ranks.

Relative reconstruction errors are **0.61273** at rank 8, **0.59949** at rank 16 and **0.59232** at rank 24. The plot shows how reconstruction changes with factor count on this split. Because tuning and gene selection overlap the split, use it as a projection example rather than an independent performance estimate.

The x-axis is factor count and the y-axis is relative held-out reconstruction error. Lower values indicate closer reconstruction of the held-out expression.

## Save and review the primary outputs

Primary results use `fitted_subset.h5ad` and `activity_scores.parquet`, with source `cellType` metadata retained. Review `sampling_audit.csv`, `sample_metadata.csv`, and results in `metabolic_walkthrough/`, `factor_interpretation_k4/` and `signaling/`. Preserve cell identifiers and ordering when reusing these outputs.

## View the generated plots

The code saves PNGs using a non-interactive plotting backend. Open the files in your output directory, or display one in Jupyter after the analysis completes:

```python
import os
from pathlib import Path
from IPython.display import Image, display
result_dir = Path(os.environ.get("SCKNIFE_TUTORIAL_OUTPUT", "tutorial/outputs_10000"))
display(Image(filename=str(result_dir / "rank_sensitivity.png")))
```

[Tutorial overview](../tutorial-10000.md) · [Previous: Cell–cell signaling](10000-signaling.md)
