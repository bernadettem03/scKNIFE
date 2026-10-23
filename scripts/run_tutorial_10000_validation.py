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


