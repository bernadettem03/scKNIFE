"""Descriptive signaling for the fitted 1,000 cells; no patient or condition tests."""
import os,time,json,hashlib,resource,textwrap
from pathlib import Path
os.environ['MPLBACKEND']='Agg'
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scknife.signaling import activity_lr as lr
from scknife.signaling import common
start=time.perf_counter()
root=Path.cwd();base=root/'tutorial/outputs_focused';out=base/'signaling';out.mkdir(exist_ok=True)
db_path=Path(os.environ['SCKNIFE_INTERACTION_DB'])
for mod in (lr,common):
 mod.ACTIVITY_PATH=base/'activity_scores.parquet';mod.CELLCHAT_DB_PATH=db_path
adata=sc.read_h5ad(base/'fitted_subset.h5ad')
assert adata.n_obs==1000
metadata=adata.obs.copy();key='bgnmf_cell_type_highest_as'
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
pd.DataFrame({'candidate_cell_type':groups,'cells':counts}).to_csv(out/'group_counts.csv',index=False)
pairs.to_csv(out/'matched_interaction_components.csv',index=False)
fig,ax=plt.subplots(figsize=(7,6));im=ax.imshow(comm,cmap='viridis',aspect='auto')
ax.set_xticks(range(len(groups)));ax.set_xticklabels([textwrap.fill(g,16) for g in groups],rotation=35,ha='right')
ax.set_yticks(range(len(groups)));ax.set_yticklabels([textwrap.fill(g,20) for g in groups])
ax.set_xlabel('Receiver candidate cell type');ax.set_ylabel('Sender candidate cell type')
ax.set_title('Descriptive communication scores · 1,000 cells');fig.colorbar(im,ax=ax,label='Weighted activity score')
fig.tight_layout();fig.savefig(out/'signaling_heatmap.png',dpi=160,bbox_inches='tight');plt.close(fig)
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
ax.set_xticks(range(len(contexts)));ax.set_xticklabels([textwrap.fill(a,17)+'\n→\n'+textwrap.fill(b,17) for a,b in contexts.index],fontsize=8)
ax.set_yticks(range(len(totals)));ax.set_yticklabels([a+' → '+b for a,b in totals.index],fontsize=9);ax.invert_yaxis()
ax.set_xlabel('Sender → receiver candidate cell types');ax.set_ylabel('Ligand → receptor component')
ax.set_title('Top 10 ligand–receptor components in 6 strongest group pairs')
fig.colorbar(scat,ax=ax,label='log(1 + weighted score)')
fig.text(.5,.01,'Bubble area is proportional to weighted score; no significance test is shown.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.04,1,1));fig.savefig(out/'signaling_lr_bubbles.png',dpi=160,bbox_inches='tight');plt.close(fig)
report={'cells':1000,'groups':len(groups),'group_counts':dict(zip(groups,map(int,counts))),'interaction_rows_after_annotation_filter':len(db),'matched_database_rows':int(pairs.db_row.nunique()),'matched_lr_components':len(pairs),'coverage':coverage,'vectorization_max_abs_error':error,'detail_totals_match_matrix':True,'mode':'descriptive pooled activity scoring; no permutations, FDR, patient or condition comparisons','patient_ids_available':'patient_id' in metadata,'interaction_file_sha256':hashlib.sha256(db_path.read_bytes()).hexdigest(),'elapsed_seconds':time.perf_counter()-start,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'versions':common.software_versions()}
(out/'signaling_check_report.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:v for k,v in report.items() if k not in ['coverage','versions']},indent=2))
