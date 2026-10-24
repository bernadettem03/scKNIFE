# Cell–cell signaling

[Tutorial overview](../tutorial-10000.md) · [Previous: Factor interpretation](10000-factors.md) · [Optional: rank sensitivity and projection](10000-validation.md)

Explore activity-based signaling in seven steps: **communication overview → sender/receiver roles → ligand–receptor contributions → permutation FDR → network organization → VEGF → On-minus-Pre differences**.

Run the blocks in order from the repository root after fitting and activity reconstruction. All analyses use the same 10,000 cells and source `cellType` labels. Signaling scores summarize reconstructed activities using directional ligand–receptor definitions; they are not CellChat probabilities. [Download the complete runnable script](../downloads/signaling_walkthrough.py).

Set `SCKNIFE_INTERACTION_DB` to the directional interaction export described in [Input files](../tutorial-inputs.md). The interaction export supplies directional ligand–receptor definitions alongside the fitted activities.

## Load and validate signaling inputs

```python
"""Descriptive signaling for the fitted 10,000 cells; no patient or condition tests."""
import os,time,json,hashlib,resource,textwrap
from pathlib import Path
os.environ['MPLBACKEND']='Agg'
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# Consistent, readable figure typography.
plt.rcParams.update({"font.size": 12, "axes.titlesize": 14,
                     "axes.labelsize": 12, "xtick.labelsize": 11,
                     "ytick.labelsize": 11, "legend.fontsize": 11,
                     "figure.titlesize": 15, "axes.spines.top": False,
                     "axes.spines.right": False})
from scknife.signaling import activity_lr as lr
from scknife.signaling import common
start=time.perf_counter()
root=Path.cwd();base=Path(os.environ.get('SCKNIFE_TUTORIAL_OUTPUT', root/'tutorial/outputs_10000'));out=base/'signaling';out.mkdir(exist_ok=True)
db_path=Path(os.environ['SCKNIFE_INTERACTION_DB'])
for mod in (lr,common):
 mod.ACTIVITY_PATH=base/'activity_scores.parquet';mod.CELLCHAT_DB_PATH=db_path
adata=sc.read_h5ad(base/'fitted_subset.h5ad')
assert adata.n_obs==10000
metadata=adata.obs.copy();key='cellType'
assert key in metadata and not metadata[key].isna().any()
db,wanted=lr.load_interactions()
_,nodes,cells=lr.activity_schema()
assert cells==list(adata.obs_names)
values,genes,coverage=lr.stream_activity_subset(wanted,nodes,cells,256)
groups=sorted(metadata[key].astype(str).unique());assert len(groups)>1
labels=metadata[key].astype(str).map({g:i for i,g in enumerate(groups)}).to_numpy(dtype=np.int16)
means,counts=lr.group_means(values,labels,len(groups));specificity=means/(means.mean(axis=1)[:,None]+1e-10)
pairs=lr.build_pairs(db,genes,means);error=lr.validate_vectorization(means,specificity,pairs)
comm=lr.communication_matrix(means,specificity,pairs);assert np.isfinite(comm).all() and (comm>=0).all()
details=lr.interaction_details(means,specificity,pairs,groups)
reconstructed=details.pivot_table(index='source',columns='target',values='weighted_score',aggfunc='sum',fill_value=0).reindex(index=groups,columns=groups,fill_value=0)
assert np.allclose(comm,reconstructed.to_numpy(),rtol=1e-10,atol=1e-10)
pd.DataFrame(comm,index=groups,columns=groups).to_csv(out/'communication_matrix.csv')
details.to_csv(out/'interaction_details.csv.gz',index=False,compression='gzip')
pd.DataFrame({'reference_cell_type':groups,'cells':counts}).to_csv(out/'group_counts.csv',index=False)
pairs.to_csv(out/'matched_interaction_components.csv',index=False)
```

The code aligns cells with the saved fit and checks that interaction scores sum to the group matrix. Scores combine group means with group-normalized specificity. Check the group counts when interpreting small populations.

## 1. Sender–receiver activity heatmap

```python
fig,ax=plt.subplots(figsize=(7,6));im=ax.imshow(comm,cmap='viridis',aspect='auto')
ax.set_xticks(range(len(groups)));ax.set_xticklabels([textwrap.fill(g,16) for g in groups],rotation=35,ha='right')
ax.set_yticks(range(len(groups)));ax.set_yticklabels([textwrap.fill(g,20) for g in groups])
ax.set_xlabel('Receiver source cell type');ax.set_ylabel('Sender source cell type')
ax.set_title('Descriptive communication scores · 10,000 cells');fig.colorbar(im,ax=ax,label='Weighted activity score')
fig.tight_layout();fig.savefig(out/'signaling_heatmap.png',dpi=160,bbox_inches='tight');plt.close(fig)
```

```{image} ../assets/tutorial-10000/signaling-seven/signaling_heatmap.png
:class: compact-scientific-plot
:align: center
:alt: Sender–receiver activity scores
```

Rows are sender cell types and columns are receivers. Color shows summed interaction scores, including self-links. Use the FDR plot in step 4 to inspect statistical support.

## 2. Outgoing and incoming signaling strength

```python
from scknife.signaling import plotting as plots

# Consistent, readable figure typography.
plt.rcParams.update({"font.size": 12, "axes.titlesize": 14,
                     "axes.labelsize": 12, "xtick.labelsize": 11,
                     "ytick.labelsize": 11, "legend.fontsize": 11,
                     "figure.titlesize": 15, "axes.spines.top": False,
                     "axes.spines.right": False})
plots.NODE_LABEL_SIZE_PT = 11
plots.TOP_LR_PAIRS = 15
plots.OUTPUT_DPI = 180
comm = pd.DataFrame(comm, index=groups, columns=groups)
general = out / 'general'

plots.plot_sender_receiver(comm, general)
# plots.plot_flow(comm, general)  # Optional dominant-flow view.
```

```{image} ../assets/tutorial-10000/signaling-seven/01_sender_receiver.png
:class: compact-scientific-plot
:align: center
:alt: Outgoing and incoming signaling strength
```

Outgoing strength sums each row of the communication matrix; incoming strength sums each column. Both include self-links. These summarize scores rather than cell counts and identify prominent senders and receivers under the scoring model.

## 3. Ligand–receptor contributions

```python
# Aggregate database components with identical ligand/receptor labels for display.
totals=details.groupby(['ligand','receptor'],sort=True).weighted_score.sum().sort_values(ascending=False,kind='stable').head(10)
contexts=details.groupby(['source','target'],sort=True).weighted_score.sum().sort_values(ascending=False,kind='stable').head(6)
agg=details.groupby(['ligand','receptor','source','target'],sort=True).weighted_score.sum()
rows=[]
for y,(lig,rec) in enumerate(totals.index):
 for x,(sender,receiver) in enumerate(contexts.index):
  score=float(agg.get((lig,rec,sender,receiver),0))
  rows.append(dict(ligand=lig,receptor=rec,source=sender,target=receiver,score=score,x=x,y=y))
bubbles=pd.DataFrame(rows);bubbles.to_csv(out/'selected_lr_bubbles.csv',index=False)
fig,ax=plt.subplots(figsize=(10,8));positive=bubbles[bubbles.score>0];sizes=280*positive.score/positive.score.max()
scat=ax.scatter(positive.x,positive.y,s=sizes,c=np.log1p(positive.score),cmap='viridis',edgecolors='gray',linewidths=.3)
ax.set_xticks(range(len(contexts)));ax.set_xticklabels([textwrap.fill(a,17)+'\n→\n'+textwrap.fill(b,17) for a,b in contexts.index],fontsize=11)
ax.set_yticks(range(len(totals)));ax.set_yticklabels([a+' → '+b for a,b in totals.index],fontsize=11);ax.invert_yaxis()
ax.set_xlabel('Sender → receiver source cell types');ax.set_ylabel('Ligand → receptor component')
ax.set_title('Top 10 ligand–receptor components in 6 strongest group pairs')
fig.colorbar(scat,ax=ax,label='log(1 + weighted score)')
fig.text(.5,.01,'Bubble area is proportional to weighted score; no significance test is shown.',ha='center',fontsize=11)
fig.tight_layout(rect=(0,.04,1,1));fig.savefig(out/'signaling_lr_bubbles.png',dpi=160,bbox_inches='tight');plt.close(fig)
# plots.plot_bubble(details, groups, general)  # Alternative group-pair overview.
# plots.plot_top_lr(details, general)  # Alternative pooled ligand–receptor ranking.
```

```{image} ../assets/tutorial-10000/signaling-seven/signaling_lr_bubbles.png
:class: compact-scientific-plot
:align: center
:alt: Ligand–receptor scores in selected sender–receiver contexts
```

The ten highest-scoring ligand–receptor components and six strongest sender–receiver contexts are selected by pooled score. Bubble area represents score and color represents log(1 + score); neither represents FDR. The plot identifies interactions contributing to the group-level summaries.

## 4. Cell-label permutation FDR

Run or resume 1,000 global cell-label permutations with observed specificity held fixed. The null preserves group sizes. Apply Benjamini–Hochberg correction across all 64 sender–receiver pairs.

```python
from pathlib import Path
import os,json,hashlib
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from scknife.signaling import activity_lr as lr, common
base=Path(os.environ.get('SCKNIFE_TUTORIAL_OUTPUT','tutorial/outputs_10000'))
out=base/'signaling'/'inference';out.mkdir(parents=True,exist_ok=True)
for mod in (lr,common):
 mod.ACTIVITY_PATH=base/'activity_scores.parquet'
 mod.CELLCHAT_DB_PATH=Path(os.environ['SCKNIFE_INTERACTION_DB'])
 mod.DATA_DIR=out
adata=sc.read_h5ad(base/'fitted_subset.h5ad');meta=adata.obs.copy()
assert len(meta)==10000 and meta.index.is_unique
print(meta[['timepoint','expansion']].value_counts(dropna=False),flush=True)
db,wanted=lr.load_interactions();_,nodes,cells=lr.activity_schema()
assert cells==list(meta.index)
values,genes,coverage=lr.stream_activity_subset(wanted,nodes,cells,256)
groups=sorted(meta.cellType.astype(str).unique());labels=meta.cellType.astype(str).map({g:i for i,g in enumerate(groups)}).to_numpy(dtype=np.int16)
means,counts=lr.group_means(values,labels,len(groups));specificity=means/(means.mean(axis=1)[:,None]+1e-10)
pairs=lr.build_pairs(db,genes,means);observed=lr.communication_matrix(means,specificity,pairs)
saved=pd.read_csv(base/'signaling/communication_matrix.csv',index_col=0).loc[groups,groups]
assert np.allclose(observed,saved.to_numpy())
# Cell-label null; keep observed specificity fixed, matching the package method.
sig=hashlib.sha256(labels.tobytes()+json.dumps(groups).encode()).hexdigest()
pvalues,completed=lr.run_permutations(values,labels,observed,specificity,pairs,len(groups),1000,50,False,sig)
fdr=multipletests(pvalues.ravel(),method='fdr_bh')[1].reshape(pvalues.shape)
for name,v in [('communication_pvalues',pvalues),('communication_fdr',fdr)]:pd.DataFrame(v,index=groups,columns=groups).to_csv(out/(name+'.csv'))
def heatmap(v,title,filename,label,diverging=False):
 fig,ax=plt.subplots(figsize=(9,7));cmap=plt.get_cmap('RdBu_r' if diverging else 'viridis').copy();cmap.set_bad('#dddddd')
 bound=max(float(np.nanmax(np.abs(v))),1e-12)
 im=ax.imshow(v,cmap=cmap,aspect='auto',vmin=-bound if diverging else 0,vmax=bound)
 ax.set_xticks(range(len(groups)),groups,rotation=40,ha='right');ax.set_yticks(range(len(groups)),groups)
 ax.set_xlabel('Receiver cell type');ax.set_ylabel('Sender cell type');ax.set_title(title)
 fig.colorbar(im,ax=ax,label=label);fig.tight_layout();fig.savefig(out/filename,dpi=170);plt.close(fig)


heatmap(-np.log10(fdr),'Cell-label permutation FDR · 1,000 permutations','permutation_fdr.png','−log10(BH-adjusted p-value)')
```

```{image} ../assets/tutorial-10000/signaling-seven/permutation_fdr.png
:class: compact-scientific-plot
:align: center
:alt: Cell-label permutation FDR for sender–receiver pairs
```

Color shows −log10 of the BH-adjusted p-value for the 64 sender–receiver totals. These global cell-label permutations preserve group sizes but do not block by patient. The FDR applies to these totals, not individual interactions or treatment contrasts; 1,000 permutations give a minimum unadjusted p-value of 1/1,001.

## 5. Network organization and centrality

```python
plots.plot_network(comm, general, title='Cell-type communication network')
```

```{image} ../assets/tutorial-10000/signaling-seven/06_hierarchical_network.png
:class: compact-scientific-plot
:align: center
:alt: Cell-type signaling network
```

Arrows show direction for up to 40 positive links; widths use the 2nd–98th score percentiles. Circular positions arrange the network for display. Other positive links remain in the exported matrix.

```python
# Network centrality: omit self-links; stronger edges have shorter inverse-score distances.
g=nx.DiGraph();g.add_nodes_from(groups)
for i,a in enumerate(groups):
 for j,b in enumerate(groups):
  if i!=j and observed[i,j]>0:g.add_edge(a,b,weight=float(observed[i,j]),distance=1/float(observed[i,j]))
between=nx.betweenness_centrality(g,weight='distance',normalized=True)
rank=nx.pagerank(g,weight='weight')
central=pd.DataFrame({'cell_type':groups,'outgoing_strength':[g.out_degree(a,weight='weight') for a in groups],'incoming_strength':[g.in_degree(a,weight='weight') for a in groups],'betweenness':[between[a] for a in groups],'pagerank':[rank[a] for a in groups]}).set_index('cell_type')
central.to_csv(out/'network_centrality.csv')
fig,axes=plt.subplots(2,2,figsize=(14,10))
for ax,col in zip(axes.flat,central):
 ax.barh(groups,central[col]);ax.invert_yaxis();ax.set_title(col.replace('_',' ').title());ax.set_xlabel('Score')
fig.suptitle('Descriptive centrality · self-links excluded');fig.tight_layout();fig.savefig(out/'network_centrality.png',dpi=170);plt.close(fig)
```

```{image} ../assets/tutorial-10000/signaling-seven/network_centrality.png
:class: compact-scientific-plot
:align: center
:alt: Network centrality excluding self-links
```

Centrality uses the full positive-score network, excluding self-links. Strength sums edge scores, betweenness uses inverse score as distance, and PageRank uses score as weight. Step 2 includes self-links, so its totals differ.

## 6. Pathway example: VEGF

```python
SELECTED_PATHWAY = 'VEGF'
# SELECTED_PATHWAY = 'CXCL'
# SELECTED_PATHWAY = 'TGFb'
# SELECTED_PATHWAY = 'WNT'
# SELECTED_PATHWAY = 'NOTCH'
# SELECTED_PATHWAY = 'IFN'
# SELECTED_PATHWAY = 'TNF'
frame = details.loc[details.pathway_name.astype(str) == SELECTED_PATHWAY].copy()
if frame.empty:
    raise ValueError('No interactions for ' + SELECTED_PATHWAY)
path_comm = frame.pivot_table(index='source', columns='target', values='weighted_score', aggfunc='sum', fill_value=0).reindex(index=groups, columns=groups, fill_value=0)
path_dir = base / 'signaling' / 'pathways' / SELECTED_PATHWAY

plots.plot_heatmaps(path_comm, None, path_dir, title=f'{SELECTED_PATHWAY} communication')
```

```{image} ../assets/tutorial-10000/signaling-seven/VEGF_heatmap.png
:class: compact-scientific-plot
:align: center
:alt: VEGF sender–receiver activity scores
```

This heatmap summarizes VEGF records from the interaction table using a pathway-specific color range. Compare pairs within this plot. Step 4’s FDR covers global totals, not this pathway; choose another available pathway using the commented alternatives.

```python
plots.plot_top_lr(frame, path_dir)
# plots.plot_sender_receiver(path_comm, path_dir)
# plots.plot_flow(path_comm, path_dir)
# plots.plot_bubble(frame, groups, path_dir)
# plots.plot_network(path_comm, path_dir, title=f"{SELECTED_PATHWAY} cell-type network")
```

```{image} ../assets/tutorial-10000/signaling-seven/VEGF_top_lr.png
:class: compact-scientific-plot
:align: center
:alt: Top VEGF ligand–receptor pairs
```

The ranking shows up to 15 VEGF ligand–receptor pairs by weighted score summed across cell-type pairs. It explains which interactions contribute to this pathway summary. Uncomment an alternative pathway or plotting call to explore another view.

## 7. On-minus-Pre signaling differences

Use exact timepoint labels and require at least ten cells per cell type in each condition. Pooled specificity and interaction definitions remain fixed across conditions.

```python
contrast_report = {}
def draw_contrast(column, left, right, stem, title):
 matrices=[];group_counts=[]
 for condition in [left,right]:
  mask=meta[column].astype(str).eq(condition).to_numpy()
  m,c=lr.group_means(values[mask],labels[mask],len(groups));matrices.append(lr.communication_matrix(m,specificity,pairs));group_counts.append(c)
 # Require at least ten cells in each group/condition; unsupported pairs remain missing.
 eligible=(group_counts[0]>=10)&(group_counts[1]>=10)
 diff=matrices[0]-matrices[1];diff[~(eligible[:,None]&eligible[None,:])]=np.nan
 pd.DataFrame(diff,index=groups,columns=groups).to_csv(out/(stem+'.csv'))
 pd.DataFrame({left:group_counts[0],right:group_counts[1]},index=groups).to_csv(out/(stem+'_counts.csv'))
 heatmap(diff,title,stem+'.png',left+' − '+right+' weighted score',True)
 contrast_report[stem]={'excluded_cells':int((~meta[column].astype(str).isin([left,right])).sum()),'eligible_groups':list(np.array(groups)[eligible]),'counts':{left:list(map(int,group_counts[0])),right:list(map(int,group_counts[1]))}}

draw_contrast('timepoint', 'On', 'Pre', 'on_minus_pre', 'On minus Pre · pooled activity difference')
# draw_contrast('expansion', 'E', 'NE', 'expanded_minus_nonexpanded',
#               'E minus NE · clonotype expansion groups')
```

```{image} ../assets/tutorial-10000/signaling-seven/on_minus_pre.png
:class: compact-scientific-plot
:align: center
:alt: Descriptive On-minus-Pre signaling differences
```

Positive values indicate higher pooled scores On-treatment; negative values indicate higher scores Pre-treatment. Grey pairs have fewer than ten cells in one or both conditions. This descriptive contrast has no patient adjustment or contrast p-values. The E/NE alternative groups clonotype expansion, not clinical response.

## Outputs and next step

Scores and interaction tables are saved in `signaling/`; network plots in `signaling/general/`; permutation results, centrality and contrasts in `signaling/inference/`; pathway plots in `signaling/pathways/VEGF/`. For an optional extension, see [rank sensitivity and projection](10000-validation.md).
