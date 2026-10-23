from pathlib import Path
import os, json, hashlib, textwrap
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path(os.environ.get('SCKNIFE_TUTORIAL_OUTPUT','tutorial/outputs_10000'))
OUT=BASE/'metabolic_walkthrough';OUT.mkdir(exist_ok=True)
fit=sc.read_h5ad(BASE/'fitted_subset.h5ad')
activity=pd.read_parquet(BASE/'activity_scores.parquet')
assert fit.n_obs==10000 and fit.obs_names.is_unique
assert activity.columns.equals(fit.obs_names) and activity.index.is_unique
assert 'X_umap' in fit.obsm and fit.obsm['X_umap'].shape==(10000,2)
meta=fit.obs[['cellType','patient_id','timepoint']].astype(str)
G=nx.read_graphml('graphs/version_1/unified_graph.graphml')
tasks=sorted(n for n,d in G.nodes(data=True) if d.get('kind')=='task' and n in activity.index)
assert len(tasks)>0
scores=activity.loc[tasks].T
assert np.isfinite(scores.to_numpy()).all() and scores.min().min()>=0
plt.rcParams.update({'font.size':12,'axes.titlesize':13,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10,'axes.spines.top':False,'axes.spines.right':False})
def save(name):
 plt.gcf().savefig(OUT/(name+'.png'),dpi=170,bbox_inches='tight');plt.close()
def heat(table,title,name,label,cmap='viridis',signed=False):
 fig,ax=plt.subplots(figsize=(10,max(5,.48*len(table)+1)))
 limit=max(float(np.abs(table.to_numpy()).max()),1e-12)
 im=ax.imshow(table,aspect='auto',cmap=cmap,vmin=-limit if signed else 0,vmax=limit)
 ax.set_xticks(range(table.shape[1]),table.columns,rotation=35,ha='right')
 ax.set_yticks(range(len(table)),[textwrap.fill(str(x),48) for x in table.index])
 ax.set_title(title);fig.colorbar(im,ax=ax,label=label);fig.tight_layout();save(name)


means=scores.groupby(meta.cellType,observed=True).mean().T
ranked=scores.mean().sort_values(ascending=False,kind='stable')
chosen=ranked.head(8).index
# chosen=ranked.head(12).index  # Show more tasks if needed.
means.to_csv(OUT/'task_means.csv')
heat(means.loc[chosen],'Selected metabolic-task activity','task_overview','Mean reconstructed task activity')


example_tasks=[
 'ATP generation from glucose (hypoxic conditions) - glycolysis',
 'ATP regeneration from glucose (normoxic conditions) - glycolysis + krebs cycle',
 'Acetoacetate synthesis',
]
# example_tasks=ranked.head(3).index.tolist()  # Data-ranked alternative.
assert set(example_tasks).issubset(scores.columns)
xy=fit.obsm['X_umap']
fig,axes=plt.subplots(1,3,figsize=(14,4.5))
for ax,task in zip(axes,example_tasks):
 values=np.log1p(scores[task].to_numpy())
 lo,hi=np.quantile(values,[.02,.98])
 order=np.argsort(values,kind='stable')
 im=ax.scatter(xy[order,0],xy[order,1],c=values[order],s=3,cmap='viridis',vmin=lo,vmax=hi,rasterized=True)
 ax.set_title(textwrap.fill(task,30),fontsize=11);ax.set_xticks([]);ax.set_yticks([])
 ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
 fig.colorbar(im,ax=ax,label='log(1 + task activity)',shrink=.7)
fig.tight_layout();save('task_umaps')


EXAMPLE_TASK=example_tasks[0]
# EXAMPLE_TASK=example_tasks[1]  # Oxidative ATP regeneration.
# EXAMPLE_TASK=example_tasks[2]  # Acetoacetate synthesis.
types=sorted(meta.cellType.unique())
values=np.log1p(scores[EXAMPLE_TASK])
fig,ax=plt.subplots(figsize=(10,5))
ax.violinplot([values.loc[meta.cellType.eq(t)].to_numpy() for t in types],showmedians=True)
ax.set_xticks(range(1,len(types)+1),[f'{t}\n(n={meta.cellType.eq(t).sum():,})' for t in types],rotation=35,ha='right')
ax.set_ylabel('log(1 + reconstructed task activity)')
ax.set_title(textwrap.fill(EXAMPLE_TASK,65));fig.tight_layout();save('task_distribution')


# Traverse only existing task–reaction and reaction–gene edges.
reactions=sorted(n for n in G.neighbors(EXAMPLE_TASK) if G.nodes[n].get('kind')=='reaction' and n in activity.index)
assert reactions, 'No linked reactions available for this task'
genes=sorted({n for r in reactions for n in G.neighbors(r) if G.nodes[n].get('kind')=='gene' and n in activity.index})
assert genes, 'No linked genes available for this task'
def top_nodes(names,n):
 return activity.loc[names].mean(axis=1).sort_values(ascending=False,kind='stable').head(n).index.tolist()
shown_reactions=top_nodes(reactions,6)
shown_genes=top_nodes(genes,8)
# shown_reactions=top_nodes(reactions,10)
# shown_genes=top_nodes(genes,12)
linked=[EXAMPLE_TASK]+shown_reactions+shown_genes
linked_means=activity.loc[linked].T.groupby(meta.cellType,observed=True).mean().T
logged=np.log1p(linked_means)
relative=logged.sub(logged.mean(axis=1),axis=0).div(logged.std(axis=1,ddof=0).clip(lower=1e-12),axis=0)
relative.index=[f'{G.nodes[n].get("kind")}: {n}' for n in linked]
heat(relative,'Task and connected reaction/gene activity','linked_support','Row z-score of log(1 + mean activity)','RdBu_r',True)
linked_means.to_csv(OUT/'linked_node_means.csv')
edges=[{'source':EXAMPLE_TASK,'target':r,'relation':'task–reaction'} for r in reactions]
edges += [{'source':r,'target':g,'relation':'reaction–gene'} for r in reactions for g in G.neighbors(r) if g in genes]
pd.DataFrame(edges).to_csv(OUT/'task_reaction_gene_links.csv',index=False)


CELL_TYPE='Cancer_cell'
# CELL_TYPE='T_cell'
# CELL_TYPE='Myeloid_cell'
mask=meta.cellType.eq(CELL_TYPE) & meta.timepoint.isin(['Pre','On'])
counts=meta.loc[mask].groupby('timepoint').size()
assert all(counts.get(t,0)>=10 for t in ['Pre','On'])
condition_means=scores.loc[mask].groupby(meta.loc[mask,'timepoint']).mean()
delta=condition_means.loc['On']-condition_means.loc['Pre']
selected_delta=delta.loc[delta.abs().sort_values(ascending=False,kind='stable').head(8).index].sort_values()
fig,ax=plt.subplots(figsize=(10,6))
ax.barh(range(len(selected_delta)),selected_delta,color=np.where(selected_delta>=0,'#b85b63','#367fa5'))
ax.set_yticks(range(len(selected_delta)),[textwrap.fill(t,48) for t in selected_delta.index])
ax.axvline(0,color='#333333',lw=.8);ax.set_xlabel('Mean reconstructed activity: On − Pre')
ax.set_title(f'{CELL_TYPE} · descriptive task differences\nPre n={counts["Pre"]:,}; On n={counts["On"]:,}')
fig.tight_layout();save('ranked_differences')
condition_means.T.assign(On_minus_Pre=delta).to_csv(OUT/'condition_task_means.csv')
meta.loc[mask].groupby(['patient_id','timepoint']).size().unstack(fill_value=0).to_csv(OUT/'condition_patient_counts.csv')
summary={'cells':fit.n_obs,'tasks':len(tasks),'example_tasks':example_tasks,'example_task':EXAMPLE_TASK,'cell_type':CELL_TYPE,'linked_reactions':len(reactions),'linked_genes':len(genes),'condition_cells':counts.to_dict(),'cell_id_sha256':hashlib.sha256('\n'.join(fit.obs_names).encode()).hexdigest()}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2));print(summary)
