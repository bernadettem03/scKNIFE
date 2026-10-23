# Robustness checks and exported results

[Tutorial overview](../tutorial-10000.md) · [Previous: Cell–cell signaling](10000-signaling.md)

**Prerequisite:** [model fitting](10000-setup.md) for projection; complete the topic pages before reviewing their exported results. Run the command or its Python block, then review the exported results.

```bash
python scripts/run_tutorial_10000_validation.py
```

Split the selected sample into **8,000 training and 2,000 held-out cells**, seed 42. Fit ranks 8, 16 and 24 on training cells and project held-out expression onto each fixed `Rg @ W` basis using nonnegative least squares. The primary biological analysis remains the separate 16-factor fit on all 10,000 cells.

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
rows=[]; projections=[]
for rank in (8,16,24):
    model=run_full_pipeline(train.copy(),G,k=rank,lambda_reg=1,beta_W=.01,beta_H=.01,n_iter=60,tol=0,random_state=42,verbose=False)
    basis=np.asarray(model.uns['Rg'] @ model.uns['W_graph'])
    genes=model.var_names
    expression=test[:,genes].X
    expression=expression.toarray() if hasattr(expression,'toarray') else np.asarray(expression)
    predicted=np.column_stack([nnls(basis,x)[0] for x in expression])
    reconstructed=(basis@predicted).T
    error=np.linalg.norm(expression-reconstructed)/np.linalg.norm(expression)
    rows.append({'rank':rank,'training_cells':8000,'heldout_cells':2000,'relative_reconstruction_error':float(error)})
    if rank==16:
        projections=pd.DataFrame(predicted.T,index=test.obs_names,columns=[f'H_projected_{i}' for i in range(rank)])
    plt.close('all')
pd.DataFrame(rows).to_csv(OUT/'rank_sensitivity.csv',index=False)
projections.to_csv(OUT/'heldout_projected_factors.csv')
fig,ax=plt.subplots(figsize=(6,4));ax.plot([r['rank'] for r in rows],[r['relative_reconstruction_error'] for r in rows],marker='o')
ax.set_xticks([8,16,24]);ax.set_xlabel('Number of factors');ax.set_ylabel('Held-out relative reconstruction error')
ax.set_title('Exploratory rank sensitivity · 2,000 held-out cells')
fig.tight_layout();fig.savefig(OUT/'rank_sensitivity.png',dpi=160);plt.close(fig)
report={'elapsed_seconds':time.perf_counter()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'training_cells':8000,'heldout_cells':2000,'genes_selected_on_full_tutorial_sample':True,'purpose':'exploratory reconstruction diagnostic, not unbiased generalization or biological rank selection'}
(OUT/'validation_report.json').write_text(json.dumps(report,indent=2))
print(rows,report)

```

```{image} ../assets/tutorial-10000/rank_sensitivity.png
:class: compact-scientific-plot
:align: center
:alt: Relative held-out reconstruction error across the three tested ranks.
```

Relative held-out reconstruction error across the three tested ranks.

Errors decrease from **0.61294** at rank 8 to **0.59972** at rank 16 and **0.59314** at rank 24. Because gene selection used the full sample, this is an exploratory reconstruction diagnostic, not an unbiased generalization benchmark or proof of an optimal biological rank.

The diagnostic plot has factor count on the x-axis and relative held-out reconstruction error on the y-axis; lower error means better reconstruction for this split. It does not establish better cell identities.

## Save and review the primary outputs

Primary results use `fitted_subset.h5ad` and `activity_scores.parquet`, with source `cellType` metadata retained. Review `sampling_audit.csv`, `sample_metadata.csv`, and results in `metabolic_walkthrough/`, `factor_interpretation_k4/` and `signaling/`. Preserve cell identifiers and ordering when reusing these outputs.

## Clean-run verification

The complete command sequence was tested in a fresh Python environment and isolated output directory, starting from conversion of the original R count matrix. The directional interaction table was freshly exported from the pinned database snapshot. All six stages completed, including 1,000 fresh signaling permutations. Cell identifiers matched the documented 10,000-cell sample, all 32 displayed tutorial figures were generated, and graph checksums remained unchanged.

| Stage | Wall time (seconds) |
| --- | ---: |
| Fit and reconstruct activities | 101.63 |
| Cell-type retention | 1.85 |
| Metabolic analysis | 8.72 |
| Factor interpretation | 12.22 |
| Signaling, including 1,000 permutations | 44.52 |
| Rank/projection diagnostic | 102.02 |

Times are specific to the tested machine and include process startup and input loading. Peak process memory was approximately 3.53 GiB for fitting and 2.31 GiB for the rank/projection diagnostic; other stages were not separately memory-profiled. Installation and raw-input conversion are excluded from the table.

Download the {download}`clean-run verification <../downloads/tutorial_clean_run_verification.json>`, {download}`fit report <../downloads/tutorial_10000_run_report.json>` and {download}`projection report <../downloads/tutorial_10000_validation_report.json>`.

CI checks that the six command scripts match the documented Python code, alongside package tests, supplied-graph validation and a strict documentation build. The external-data biological run is verified separately; CI's bundled synthetic 1,000-cell fixture remains a regression check.

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
