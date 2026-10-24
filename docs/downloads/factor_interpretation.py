from pathlib import Path
import os, json, hashlib, textwrap, time
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from matplotlib.patches import Patch
START=time.time()
BASE=Path(os.environ.get('SCKNIFE_TUTORIAL_OUTPUT','tutorial/outputs_10000'))
OUT=BASE/'factor_interpretation_k4';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':12,'axes.titlesize':14,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10,'axes.spines.top':False,'axes.spines.right':False})
a=sc.read_h5ad(BASE/'fitted_subset.h5ad')
assert a.n_obs==10000 and a.obs_names.is_unique
nodes=np.asarray(a.uns['node_names']).astype(str)
W=np.asarray(a.uns['W_graph'],dtype=float)
hcols=sorted([c for c in a.obs if c.startswith('H_graph_')],key=lambda c:int(c.rsplit('_',1)[1]))
H=a.obs[hcols].to_numpy(float)
N_FACTORS = H.shape[1]
assert H.shape[0] == 10000 and W.shape == (len(nodes), N_FACTORS)
assert np.isfinite(W).all() and np.isfinite(H).all() and W.min()>=0 and H.min()>=0
meta=a.obs[['cellType','patient_id','timepoint']].astype(str)
types=sorted(meta.cellType.unique())
G=nx.read_graphml('graphs/version_1/unified_graph.graphml')
kinds=np.array([str(G.nodes[n].get('kind','')) for n in nodes])
keep=(~np.isin(kinds,['cell_type','cell_state'])) & (W.sum(axis=1)>0)
W_ELIGIBLE=pd.DataFrame(W[keep],index=nodes[keep],columns=[str(i) for i in range(N_FACTORS)])
NODE_METADATA_ELIGIBLE=pd.DataFrame({'node_name':nodes[keep],'node_kind':kinds[keep]})
def row_z(x):
 x=np.log1p(np.asarray(x,dtype=float))
 return (x-x.mean(axis=1,keepdims=True))/np.maximum(x.std(axis=1,keepdims=True),1e-12)
def save(name):
 plt.gcf().savefig(OUT/(name+'.png'),dpi=170,bbox_inches='tight');plt.close()
def heat(table,title,name,label='Row z-score of log(1 + mean activity)',signed=True):
 fig,ax=plt.subplots(figsize=(9,max(4,.33*len(table)+1.5)))
 lim=max(float(np.max(np.abs(table.to_numpy()))),1e-12)
 im=ax.imshow(table,aspect='auto',cmap='RdBu_r' if signed else 'viridis',vmin=-lim if signed else 0,vmax=lim)
 ax.set_xticks(range(table.shape[1]),table.columns,rotation=40,ha='right')
 ax.set_yticks(range(len(table)),table.index);ax.set_title(title)
 fig.colorbar(im,ax=ax,label=label);fig.tight_layout();save(name)


factor_means=pd.DataFrame(H,index=a.obs_names,columns=hcols).groupby(meta.cellType).mean().reindex(types).T
heat(factor_means,'Mean H activity by source cell type','factor_usage','Mean H activity',False)
factor_means.to_csv(OUT/'factor_usage.csv')


SEED=42;BOOTSTRAPS=200;RESAMPLE_FRACTION=.8
K_VALUES=range(2,min(13,N_FACTORS))
W_Z=row_z(W_ELIGIBLE)
rng=np.random.default_rng(SEED)
sums={k:np.zeros((N_FACTORS,N_FACTORS)) for k in K_VALUES}
for b in range(BOOTSTRAPS):
 rows=rng.choice(len(W_Z),int(.8*len(W_Z)),replace=True)
 d=np.clip(1-np.corrcoef(W_Z[rows].T),0,2);np.fill_diagonal(d,0)
 for k in K_VALUES:
  lab=AgglomerativeClustering(n_clusters=k,metric='precomputed',linkage='average').fit_predict(d)
  sums[k]+=lab[:,None]==lab[None,:]
records=[];solutions={}
for k in K_VALUES:
 c=sums[k]/BOOTSTRAPS;d=1-c;np.fill_diagonal(d,0)
 tree=linkage(squareform(d,checks=False),method='average')
 lab=fcluster(tree,t=k,criterion='maxclust')-1
 sizes=np.bincount(lab);upper=c[np.triu_indices(N_FACTORS,1)]
 records.append({'k':k,'actual_groups':len(sizes),'PAC':np.mean((upper>.1)&(upper<.9)),
  'silhouette':silhouette_score(d,lab,metric='precomputed'),'min_size':sizes.min(),'max_size':sizes.max()})
 solutions[k]=(c,lab,tree)
metrics=pd.DataFrame(records)
# Review consensus diagnostics, then use the chosen four-group resolution.
eligible=metrics[(metrics.min_size>=2)&(metrics.actual_groups==metrics.k)]
assert not eligible.empty
SELECTED_K=4  # Chosen tutorial resolution; fitted model rank is selected before the full fit.
assert SELECTED_K in eligible.k.values
metrics.to_csv(OUT/'k_diagnostics.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(10,4))
for ax,col,title in zip(axes,['PAC','silhouette'],['Ambiguous consensus fraction · lower is better','Consensus silhouette · higher is better']):
 ax.plot(metrics.k,metrics[col],'o-',color='#247b87')
 ax.scatter(eligible.k,eligible[col],s=110,facecolors='none',edgecolors='#d7832c',label='No singleton groups')
 ax.axvline(SELECTED_K,color='#9b456b',ls='--',label=f'Selected k={SELECTED_K}')
 ax.set_xlabel('Number of factor groups (k)');ax.set_title(title,fontsize=12);ax.set_xticks(list(K_VALUES))
axes[0].legend(fontsize=9);fig.tight_layout();save('choose_k')
print(metrics.to_string(index=False),flush=True);print('SELECTED_K',SELECTED_K,flush=True)


consensus,raw,tree=solutions[SELECTED_K]
raw_order=sorted(set(raw),key=lambda l:np.flatnonzero(raw==l)[0])
labels=np.array([raw_order.index(l) for l in raw])
groups=[f'CC{i+1}' for i in range(SELECTED_K)]
ngroups=[f'NG{i+1}' for i in range(SELECTED_K)]
COMPONENT_ASSIGNMENTS=pd.DataFrame({'component':[str(i) for i in range(N_FACTORS)],'consensus_cluster':labels})
COMPONENT_ASSIGNMENTS.to_csv(OUT/'factor_groups.csv',index=False)
order=np.lexsort((np.arange(N_FACTORS),labels))
ticks=[f'H{i} · {groups[labels[i]]}' for i in order]
heat(pd.DataFrame(consensus[np.ix_(order,order)],index=ticks,columns=ticks),'Bootstrap factor co-clustering','consensus','Co-clustering fraction',False)


sensitivity=[]
for k in [2,3,4]:
 lab=solutions[k][1]
 ordered_labels=sorted(set(lab),key=lambda v:np.flatnonzero(lab==v)[0])
 sensitivity.append([ordered_labels.index(v)+1 for v in lab])
fig,ax=plt.subplots(figsize=(10,3))
ax.imshow(sensitivity,cmap='tab10',aspect='auto',vmin=1,vmax=10)
for i,row in enumerate(sensitivity):
 for j,v in enumerate(row):ax.text(j,i,str(v),ha='center',va='center',color='white')
ax.set_xticks(range(N_FACTORS),[f'H{i}' for i in range(N_FACTORS)],rotation=45)
ax.set_yticks(range(3),['k=2','k=3','k=4']);ax.set_title('Factor membership at nearby resolutions')
fig.tight_layout();save('resolution_sensitivity')


DECAY_TAU=5.;DECAY_CUMULATIVE_MASS=.95;MAX_NODES_PER_CLUSTER_KIND=20;TOP_N=20
from scknife.tutorial_helpers import (
    score_nodes, assign_nodes_to_strongest_cluster,
    rank_nodes_with_exponential_decay, select_mixed_decay_budget,
)

ALL_NODE_SCORES=score_nodes(W_Z,NODE_METADATA_ELIGIBLE,COMPONENT_ASSIGNMENTS,SELECTED_K)
assigned=assign_nodes_to_strongest_cluster(ALL_NODE_SCORES)
ranked,decay_audit=rank_nodes_with_exponential_decay(assigned,DECAY_TAU,DECAY_CUMULATIVE_MASS,MAX_NODES_PER_CLUSTER_KIND)
local,local_audit=rank_nodes_with_exponential_decay(ALL_NODE_SCORES[ALL_NODE_SCORES.enrichment_delta_z>0],DECAY_TAU,DECAY_CUMULATIVE_MASS,MAX_NODES_PER_CLUSTER_KIND)
selected=select_mixed_decay_budget(ranked,local,TOP_N,SELECTED_K)
assert selected.node_name.is_unique and len(selected)==TOP_N*SELECTED_K
assert selected.enrichment_delta_z.gt(0).all()
selected.to_csv(OUT/'selected_nodes.csv',index=False);decay_audit.to_csv(OUT/'decay_audit.csv',index=False)


wmean=np.column_stack([W_ELIGIBLE.to_numpy()[:,labels==g].mean(axis=1) for g in range(SELECTED_K)])
wz=pd.DataFrame(row_z(wmean),index=W_ELIGIBLE.index,columns=groups)
rows=selected.sort_values(['consensus_cluster','cluster_rank'])
row_groups=rows.consensus_cluster.to_numpy()
row_names=[textwrap.fill(f'{ngroups[g]}: {n} [{kind}]',52)
           for g,n,kind in zip(row_groups,rows.node_name,rows.node_kind)]
def combined_node_heatmap(values,columns,title,filename,color_label,factor_columns=False):
 fig,ax=plt.subplots(figsize=(12,max(6,.36*len(rows)+2)))
 limit=max(float(np.abs(values).max()),1e-12)
 im=ax.imshow(values,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
 ax.set_yticks(range(len(rows)),row_names,fontsize=10)
 ax.set_xticks(range(len(columns)),columns,rotation=45,ha='right')
 for boundary in np.flatnonzero(np.diff(row_groups)!=0)+.5:
  ax.axhline(boundary,color='#273e50',lw=1.2)
 if factor_columns:
  for boundary in np.flatnonzero(np.diff(labels[order])!=0)+.5:
   ax.axvline(boundary,color='#273e50',lw=1.2)
  top=ax.secondary_xaxis('top')
  top.set_xticks([np.flatnonzero(labels[order]==g).mean() for g in range(SELECTED_K)],groups)
  top.tick_params(length=0)
 ax.set_title(title,pad=30 if factor_columns else 12)
 fig.colorbar(im,ax=ax,label=color_label,fraction=.025,pad=.03)
 fig.tight_layout();save(filename)
combined_node_heatmap(wz.loc[rows.node_name].to_numpy(),groups,
 'Representative biological nodes · NG1–NG4','node_dictionary_all',
 'Node-wise z-score of log(1 + mean W)')
combined_node_heatmap(W_Z[W_ELIGIBLE.index.get_indexer(rows.node_name)][:,order],
 [f'H{i}' for i in order],'Biological nodes across factors · NG1–NG4',
 'node_factors_all','Node-wise z-score of log(1 + W)',factor_columns=True)
composition=pd.crosstab(selected.consensus_cluster,selected.node_kind);composition.index=ngroups
ax=composition.plot.bar(stacked=True,figsize=(9,5),colormap='tab20',rot=0)
ax.set_ylabel('Selected nodes');ax.set_xlabel('Node group');ax.set_title('Biological node-kind composition')
ax.legend(title='Node kind',bbox_to_anchor=(1.02,1),loc='upper left',fontsize=10);plt.tight_layout();save('composition')


# Select one entity type for the worked example.
entity_types = ['geneset']
# entity_types.append('gene')
# entity_types.append('GO_term')
# entity_types.append('reaction')
# entity_types.append('task')
# entity_types.append('cellchat_entity')
for kind in entity_types:
 rows=ranked[ranked.decay_selected & ranked.node_kind.eq(kind)]
 rows=rows.sort_values(['consensus_cluster','category_rank'])
 values=W_Z[W_ELIGIBLE.index.get_indexer(rows.node_name)][:,order]
 row_groups=rows.consensus_cluster.to_numpy()
 names=[textwrap.fill(f'{ngroups[g]}: {n}',52) for g,n in zip(row_groups,rows.node_name)]
 fig,ax=plt.subplots(figsize=(12,max(6,.36*len(rows)+2)))
 limit=max(float(np.abs(values).max()),1e-12)
 im=ax.imshow(values,aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
 ax.set_yticks(range(len(rows)),names,fontsize=10)
 ax.set_xticks(range(N_FACTORS),[f'H{i}' for i in order],rotation=45,ha='right')
 for boundary in np.flatnonzero(np.diff(row_groups)!=0)+.5:
  ax.axhline(boundary,color='#273e50',lw=1.2)
 for boundary in np.flatnonzero(np.diff(labels[order])!=0)+.5:
  ax.axvline(boundary,color='#273e50',lw=1.2)
 top=ax.secondary_xaxis('top')
 centers=[np.flatnonzero(labels[order]==g).mean() for g in range(SELECTED_K)]
 top.set_xticks(centers,groups);top.tick_params(length=0)
 ax.set_title(f'{kind.replace("_"," ")} · all four programs',pad=30)
 ax.set_xlabel('Factors ordered by factor group')
 fig.colorbar(im,ax=ax,label='Node-wise z-score of log(1 + W)',fraction=.025,pad=.03)
 fig.tight_layout();save(f'category_all_{kind}')


profiles=np.vstack([W[np.isin(nodes,selected.loc[selected.consensus_cluster==g,'node_name'])].mean(axis=0) for g in range(SELECTED_K)])
node_raw=H@profiles.T
factor_raw=np.column_stack([H[:,labels==g].mean(axis=1) for g in range(SELECTED_K)])
node_mean=pd.DataFrame(node_raw,index=a.obs_names,columns=ngroups).groupby(meta.cellType).mean().reindex(types).T
factor_mean=pd.DataFrame(factor_raw,index=a.obs_names,columns=groups).groupby(meta.cellType).mean().reindex(types).T
node_type=pd.DataFrame(row_z(node_mean),index=ngroups,columns=types)
factor_type=pd.DataFrame(row_z(factor_mean),index=groups,columns=types)
heat(node_type,'Reconstructed node-group activity','node_activity')
heat(factor_type,'Shared W–H factor-group activity','factor_activity')
node_mean.to_csv(OUT/'node_group_means.csv');factor_mean.to_csv(OUT/'factor_group_means.csv')


node_factor=pd.DataFrame(np.vstack([wz.loc[selected.loc[selected.consensus_cluster==g,'node_name']].mean(axis=0) for g in range(SELECTED_K)]),index=ngroups,columns=groups)
fig,ax=plt.subplots(figsize=(10,6))
ys=[np.linspace(0,1,len(x)) for x in [ngroups,groups,types]]
edges=[]
for layer,table in enumerate([node_factor,factor_type]):
 for i,src in enumerate(table.index):
  for j,dst in enumerate(table.columns):
   v=float(table.iloc[i,j]);edges.append({'source':src,'target':dst,'relative_score':v})
   ax.plot([layer,layer+1],[ys[layer][i],ys[layer+1][j]],color='#b24b55' if v>=0 else '#357eaa',lw=.3+2*abs(v),alpha=.38)
for x,names in enumerate([ngroups,groups,types]):
 ax.scatter([x]*len(names),ys[x],s=85,color='#273e50',zorder=3)
 for y,name in zip(ys[x],names):ax.text(x+.04,y,name,va='center',fontsize=11)
for x,title in enumerate(['W node groups','Shared factor groups','Source cell types']):ax.text(x,1.08,title,weight='bold')
ax.set_xlim(-.05,2.5);ax.set_ylim(-.06,1.15);ax.axis('off')
ax.legend(handles=[Patch(color='#b24b55',label='Positive relative score'),Patch(color='#357eaa',label='Negative relative score')],loc='lower center',bbox_to_anchor=(.5,-.12),ncol=2)
plt.tight_layout();save('associations');pd.DataFrame(edges).to_csv(OUT/'all_associations.csv',index=False)


logged=np.log1p(node_raw);cell_z=(logged-logged.mean(axis=0))/np.maximum(logged.std(axis=0),1e-12)
ordered=meta.sort_values(['cellType','patient_id','timepoint'],kind='stable').index
positions=a.obs_names.get_indexer(ordered)
fig,axes=plt.subplots(4,1,figsize=(11,6),gridspec_kw={'height_ratios':[.25,.25,.25,3]},layout='constrained')
legend=[]
for ax,key in zip(axes[:3],['cellType','patient_id','timepoint']):
 cats=sorted(meta[key].unique());codes=pd.Categorical(meta.loc[ordered,key],categories=cats).codes
 ax.imshow(codes[None,:],aspect='auto',cmap='tab20',vmin=0,vmax=19);ax.set_yticks([0],[key]);ax.set_xticks([])
 for j,c in enumerate(cats):legend.append(Patch(color=plt.get_cmap('tab20')(j),label=f'{key}: {c}'))
im=axes[-1].imshow(cell_z[positions].T,aspect='auto',cmap='RdBu_r',vmin=-2.5,vmax=2.5)
axes[-1].set_yticks(range(SELECTED_K),ngroups);axes[-1].set_xlabel('All 10,000 cells · ordered by cell type, patient and timepoint')
fig.colorbar(im,ax=axes[-1],label='Per-group z-score of log(1 + activity)')
fig.legend(handles=legend,loc='outside lower center',ncol=4,fontsize=9);save('all_cells')


coverage=pd.crosstab(meta.cellType,meta.patient_id).reindex(types)
heat(coverage,'Cells available per patient and cell type','coverage','Number of cells',False)
# A declared descriptive display: first three nodes from NG1's ranked dictionary.
chosen=selected[selected.consensus_cluster==0].head(3).node_name.tolist()
fig,axes=plt.subplots(1,3,figsize=(12,5))
for ax,node in zip(axes,chosen):
 score=np.log1p(H@W[np.flatnonzero(nodes==node)[0]])
 ax.violinplot([score[meta.cellType.to_numpy()==t] for t in types],showmedians=True)
 ax.set_xticks(range(1,len(types)+1),types,rotation=60,ha='right');ax.set_title(textwrap.fill(node,28),fontsize=11)
 ax.set_ylabel('log(1 + reconstructed node activity)')
plt.tight_layout();save('constituents')
summary={'cells':a.n_obs,'factors':N_FACTORS,'selected_k':SELECTED_K,'group_sizes':np.bincount(labels).tolist(),'bootstraps':BOOTSTRAPS,'selected_nodes':len(selected),'seconds':time.time()-START,'cell_id_sha256':hashlib.sha256('\n'.join(a.obs_names).encode()).hexdigest()}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
