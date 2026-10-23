# Optional diagnostic: marker-based reassignment

**This is not a required scKNIFE annotation step.** It tests additional marker-based candidate labels for comparison. It does not replace `cellType`, and its outputs are not used by the main metabolic, factor or signaling analyses.

[Tutorial overview](../tutorial-10000.md) · [Previous: Set up and fit](10000-setup.md) · [Return to main tutorial](../tutorial-10000.md)

**Prerequisite:** complete [sampling and fitting](10000-setup.md). This stage reads the saved fit and full-gene source expression for those same cell identifiers.

```bash
python scripts/run_tutorial_10000_breast_annotation.py
python scripts/run_tutorial_10000_rare_refinement.py
```

**Choose the commands above or run both expanded Python blocks below in order.** They perform the same annotation steps. All page blocks run from the repository root.

Assign broad breast-context labels from measured marker expression within the fitted clusters, then refine rare mast-cell and pDC candidates using cell-level marker combinations. The rules use expression and cluster membership, not reference labels or clinical metadata. They remain development-stage heuristics, and this sample is not an independent validation cohort.

## Broad breast-context annotation code

```python
"""Run breast marker annotation without reference labels, then evaluate separately."""
import os,json,time,hashlib
from pathlib import Path
os.environ['MPLBACKEND']='Agg'
import numpy as np,pandas as pd,scanpy as sc,matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score,confusion_matrix
from scknife.celltype.breast import run_breast_celltype_assignment,BREAST_MARKERS
start=time.perf_counter();ROOT=Path.cwd();BASE=Path(os.environ.get('SCKNIFE_TUTORIAL_OUTPUT',ROOT/'tutorial/outputs_10000'));OUT=BASE/'breast_annotation';OUT.mkdir(exist_ok=True)
fit=sc.read_h5ad(BASE/'fitted_subset.h5ad');assert fit.n_obs==10000
# Load all source genes for these exact cells, not just the 2,000 fitting features.
source=sc.read_h5ad(os.environ['SCKNIFE_COUNTS'])
assert fit.obs_names.isin(source.obs_names).all()
expression=source[fit.obs_names].copy();del source
expression.obs=pd.DataFrame(index=fit.obs_names)  # No reference, patient, or clinical metadata enters the rule.
sc.pp.normalize_total(expression,target_sum=1e4);sc.pp.log1p(expression)
config={'min_score':0.5,'min_margin':0.15,'min_detection':0.1,'min_markers':2}
(OUT/'annotation_config.json').write_text(json.dumps({'thresholds':config,'markers':BREAST_MARKERS,'rule':'marker-gene cluster means standardized across clusters; fixed heuristic thresholds; CD3 gate added after marker-evidence review; development-set reference evaluation'},indent=2))
labels,decisions,scores,detection,coverage=run_breast_celltype_assignment(expression,fit.obs.leiden_clusters,**config)
assert labels.index.equals(fit.obs_names) and not labels.isna().any()
# Independent of metadata: permuting a dummy reference column cannot affect the result.
expression.obs['cellType']='unused'
check,*_=run_breast_celltype_assignment(expression,fit.obs.leiden_clusters,**config)
pd.testing.assert_series_equal(labels,check)
fit.obs[labels.name]=labels
for name,table in [('cluster_decisions',decisions),('lineage_scores',scores),('marker_detection',detection),('marker_coverage',coverage)]:table.to_csv(OUT/f'{name}.csv')
labels.to_csv(OUT/'cell_labels.csv')
# Evaluation starts only after assignments are finalized.
reference_map={'Cancer_cell':'Epithelial','Fibroblast':'Fibroblast','Endothelial_cell':'Endothelial','T_cell':'T_cell','B_cell':'B_cell','Myeloid_cell':'Myeloid','pDC':'pDC','Mast_cell':'Mast_cell'}
reference=fit.obs.cellType.astype(str).map(reference_map);assert reference.notna().all()
counts=pd.crosstab(labels,reference);counts.to_csv(OUT/'reference_counts.csv')
accepted=labels.ne('Uncertain')
pertype=[]
for name in sorted(reference.unique()):
 mask=reference.eq(name);pertype.append({'reference':name,'cells':int(mask.sum()),'correct':int((labels[mask]==name).sum()),'uncertain':int((labels[mask]=='Uncertain').sum()),'recall':float((labels[mask]==name).mean())})
pd.DataFrame(pertype).to_csv(OUT/'per_reference_recall.csv',index=False)
# Shared palette for broad reference and predicted labels; no embedding is recomputed.
fit.obs['reference_broad_lineage']=reference
categories=list(BREAST_MARKERS)+['Uncertain'];colors=list(plt.get_cmap('tab20').colors);palette=dict(zip(categories,colors));palette['Uncertain']='#bdbdbd'
fig,axes=plt.subplots(1,2,figsize=(17,6))
for ax,key,title in zip(axes,['reference_broad_lineage',labels.name],['Source reference (broad lineages)','Breast marker annotation']):
 sc.pl.umap(fit,color=key,palette=palette,ax=ax,show=False,frameon=False,title=title,size=7)
fig.suptitle('Same 10,000 cells and unchanged H-space UMAP');fig.tight_layout();fig.savefig(OUT/'breast_annotation_umap.png',dpi=180,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(11,7));im=ax.imshow(counts,aspect='auto',cmap='Blues')
ax.set_xticks(range(counts.shape[1]));ax.set_xticklabels(counts.columns,rotation=40,ha='right');ax.set_yticks(range(counts.shape[0]));ax.set_yticklabels(counts.index)
for i in range(counts.shape[0]):
 for j in range(counts.shape[1]):
  v=counts.iloc[i,j];ax.text(j,i,str(v),ha='center',va='center',color='white' if v>counts.to_numpy().max()/2 else 'black')
ax.set_xlabel('Source reference mapped to broad lineage');ax.set_ylabel('Predicted broad lineage');ax.set_title('Breast annotation versus reference · all 10,000 cells')
fig.colorbar(im,ax=ax,label='Cells');fig.tight_layout();fig.savefig(OUT/'breast_annotation_reference.png',dpi=180,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(12,7));im=ax.imshow(scores.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-3,vmax=3)
ax.set_xticks(range(scores.shape[1]));ax.set_xticklabels(scores.columns,rotation=40,ha='right');ax.set_yticks(range(scores.shape[0]));ax.set_yticklabels([f'{i}: {decisions.loc[i,"label"]}' for i in scores.index])
ax.set_title('Marker evidence by existing scKNIFE cluster');fig.colorbar(im,ax=ax,label='Mean clipped marker z-score across clusters');fig.tight_layout();fig.savefig(OUT/'breast_marker_evidence.png',dpi=180,bbox_inches='tight');plt.close(fig)
report={'cells':10000,'clusters':len(decisions),'reference_used_for_assignment':False,'same_embedding':True,'assigned_cells':int(accepted.sum()),'uncertain_cells':int((~accepted).sum()),'broad_lineage_agreement_all_cells':float((labels==reference).mean()),'broad_lineage_agreement_assigned_cells':float((labels[accepted]==reference[accepted]).mean()) if accepted.any() else None,'candidate_counts':labels.value_counts().to_dict(),'elapsed_seconds':time.perf_counter()-start,'scope':'Heuristic broad-lineage expression annotation; no malignancy inference, no independent validation dataset; original reference used only for evaluation.'}
(OUT/'breast_annotation_report.json').write_text(json.dumps(report,indent=2))
fit.write_h5ad(OUT/'annotated_10000.h5ad')
print(json.dumps(report,indent=2));print(pd.DataFrame(pertype).to_string(index=False));print(decisions.to_string())

```

## Rare-cell refinement and comparison code

```python
"""Rare-cell refinement and stromal evidence review; references used only afterward."""
import os,json,hashlib
from pathlib import Path
os.environ['MPLBACKEND']='Agg'
import numpy as np,pandas as pd,scanpy as sc,matplotlib.pyplot as plt
from scknife.celltype.rare_refinement import refine_rare_breast_labels
ROOT=Path.cwd();BASE=Path(os.environ.get('SCKNIFE_TUTORIAL_OUTPUT',ROOT/'tutorial/outputs_10000'));OUT=BASE/'rare_refinement';OUT.mkdir(exist_ok=True)
fit=sc.read_h5ad(BASE/'breast_annotation/annotated_10000.h5ad');assert fit.n_obs==10000
base=fit.obs['scknife_breast_cell_type'].copy()
a=sc.read_h5ad(os.environ['SCKNIFE_COUNTS']);expression=a[fit.obs_names].copy();del a
expression.obs=pd.DataFrame(index=fit.obs_names)
sc.pp.normalize_total(expression,target_sum=1e4);sc.pp.log1p(expression)
labels,evidence,coverage=refine_rare_breast_labels(expression,base,min_expression=0.5)
fit.obs[labels.name]=labels
labels.to_csv(OUT/'refined_labels.csv');evidence.to_csv(OUT/'cell_evidence.csv');coverage.to_csv(OUT/'marker_coverage.csv',index=False)
# Stability analysis describes sensitivity; it does not choose a threshold.
stability=[]
for threshold in (0.25,0.5,1.0):
    alternative,_,_=refine_rare_breast_labels(expression,base,min_expression=threshold)
    for lineage in ['Mast_cell','pDC']:
        selected=alternative.eq(lineage);primary=labels.eq(lineage)
        stability.append({'threshold':threshold,'lineage':lineage,'cells':int(selected.sum()),'jaccard_vs_primary':float((selected&primary).sum()/max(1,(selected|primary).sum()))})
pd.DataFrame(stability).to_csv(OUT/'threshold_stability.csv',index=False)
# Perivascular assessment does not relabel cells to agree with the reference.
stromal=['RGS5','CSPG4','MCAM','NOTCH3','MYH11','ACTA2','TAGLN','DCN','LUM','COL1A1','COL1A2','PDGFRA']
x=expression[:,stromal].X;x=x.toarray() if hasattr(x,'toarray') else np.asarray(x)
frame=pd.DataFrame(x,index=fit.obs_names,columns=stromal)
means=frame.groupby(fit.obs.leiden_clusters,observed=True).mean();detection=frame.gt(0).groupby(fit.obs.leiden_clusters,observed=True).mean()
means.to_csv(OUT/'stromal_cluster_means.csv');detection.to_csv(OUT/'stromal_cluster_detection.csv')
peri=base.eq('Perivascular');peri_support=(frame[['RGS5','CSPG4','MCAM','NOTCH3','MYH11']]>=0.5).sum(axis=1)
pd.DataFrame({'perivascular_cluster_label':peri,'perivascular_marker_count':peri_support}).to_csv(OUT/'stromal_cell_support.csv')
# Evaluate only after the expression-only refinement is fixed.
ref_map={'Cancer_cell':'Epithelial','Fibroblast':'Fibroblast','Endothelial_cell':'Endothelial','T_cell':'T_cell','B_cell':'B_cell','Myeloid_cell':'Myeloid','pDC':'pDC','Mast_cell':'Mast_cell'}
ref=fit.obs.cellType.astype(str).map(ref_map);assert ref.notna().all()
comparison=pd.crosstab(labels,ref);comparison.to_csv(OUT/'reference_counts.csv')
metrics=[]
for lineage in sorted(ref.unique()):
    predicted=labels.eq(lineage);truth=ref.eq(lineage);tp=int((predicted&truth).sum())
    metrics.append({'lineage':lineage,'reference_cells':int(truth.sum()),'predicted_cells':int(predicted.sum()),'true_positive':tp,'precision':tp/int(predicted.sum()) if predicted.any() else None,'recall':tp/int(truth.sum()),'baseline_true_positive':int((base.eq(lineage)&truth).sum())})
pd.DataFrame(metrics).to_csv(OUT/'reference_metrics.csv',index=False)
fit.obs['refinement_review']=np.where(evidence.changed,'Refined','Unchanged')
fig,axes=plt.subplots(1,3,figsize=(20,6))
all_labels=sorted(set(labels)|set(ref));palette=dict(zip(all_labels,plt.get_cmap('tab20').colors));palette['Uncertain']='#aaaaaa'
fit.obs['reference_broad_lineage']=ref
for ax,key,title in zip(axes,[labels.name,'reference_broad_lineage','refinement_review'],['Optional marker-based reassignment','Retained dataset annotation','Cells changed by marker gates']):
    sc.pl.umap(fit,color=key,ax=ax,show=False,frameon=False,title=title,size=8,palette=palette if key!='refinement_review' else {'Refined':'#d73027','Unchanged':'#dddddd'})
fig.suptitle('Same 10,000 cells; original clusters and UMAP retained');fig.tight_layout();fig.savefig(OUT/'rare_refinement_umap.png',dpi=160,bbox_inches='tight');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(16,max(6,len(means)*.4)))
for ax,table,title in zip(axes,[means,detection],['Mean normalized/log1p expression','Fraction of cells detecting each gene']):
    im=ax.imshow(table,aspect='auto',cmap='Blues');ax.set_xticks(range(len(stromal)));ax.set_xticklabels(stromal,rotation=60,ha='right');ax.set_yticks(range(len(means)));ax.set_yticklabels([str(c)+': '+str(base[fit.obs.leiden_clusters.eq(c)].iloc[0]) for c in means.index]);ax.set_title(title);fig.colorbar(im,ax=ax)
fig.tight_layout();fig.savefig(OUT/'stromal_marker_review.png',dpi=160,bbox_inches='tight');plt.close(fig)
report={'cells':10000,'changed_cells':int(evidence.changed.sum()),'agreement_before':float(base.eq(ref).mean()),'agreement_after':float(labels.eq(ref).mean()),'ambiguous_rare_cells':int(evidence.ambiguous_rare_evidence.sum()),'perivascular_cells_retained':int(peri.sum()),'perivascular_cells_with_two_markers_at_0_5':int((peri&(peri_support>=2)).sum()),'rare_metrics':[m for m in metrics if m['lineage'] in ['pDC','Mast_cell']],'reference_used_for_refinement':False,'primary_threshold':0.5,'validation_scope':'Development sample; thresholds are heuristics; references used for evaluation only.'}
(OUT/'refinement_report.json').write_text(json.dumps(report,indent=2));fit.write_h5ad(OUT/'annotated_10000_refined.h5ad')
print(json.dumps(report,indent=2));print(pd.DataFrame(stability).to_string(index=False))

```

The diagnostic UMAP is generated by the code below. Its title will be updated to say “optional reassignment” during remote synchronization.

In the optional diagnostic run, rare-cell refinement changes **63 candidate labels**; the original `cellType` column is retained. The result includes **6 mast-cell candidates and 55 pDC candidates**; five and 53, respectively, agree with the source reference. Overall broad-lineage agreement is **85.1%**. This measures agreement, not proof of identity. Epithelial labels do not establish malignancy.

**Diagnostic uncertainty, not the main annotation:** 1,055 cells remain Uncertain. Cluster 2 contains 995 cells whose Endothelial and Perivascular scores tie at the clipping limit; cluster 15 also fails the required score margin. This diagnostic assigns no confident Endothelial or Perivascular group. The main analysis still retains the dataset’s Endothelial_cell labels. This does not imply biological absence: the Uncertain cells include 696 source endothelial and 285 source fibroblast labels. Review this ambiguity before making lineage-specific claims.

## Group counts and reference comparison

| Group | Primary refined | Reference comparison |
| --- | --- | --- |
| B_cell | 496 | 426 |
| Epithelial | 2276 | 2303 |
| Fibroblast | 1379 | 1703 |
| Mast_cell | 6 | 7 |
| Myeloid | 1492 | 1490 |
| T_cell | 3241 | 3223 |
| Uncertain | 1055 | Absent |
| pDC | 55 | 55 |
| Endothelial | Absent | 793 |

| lineage | reference_cells | predicted_cells | true_positive | precision | recall |
| --- | --- | --- | --- | --- | --- |
| Mast_cell | 7 | 6 | 5 | 0.833 | 0.714 |
| pDC | 55 | 55 | 53 | 0.964 | 0.964 |

For broad-lineage comparison, source Cancer_cell maps to Epithelial, Endothelial_cell to Endothelial and Myeloid_cell to Myeloid. B_cell includes the plasma/B lineage. Unresolved cells are not reassigned to improve agreement.


## Marker rules and threshold sensitivity

First apply `run_breast_celltype_assignment` to measured expression from all source genes for these exact cells, normalized to 10,000 counts and log1p. Cluster-level marker scores standardize each gene across clusters, clip z-scores to −3 to 3, and average marker panels. Require at least two markers detected in at least 10% of cluster cells, a score of at least 0.5 and a margin of at least 0.15 over the next eligible lineage. A CD3 safeguard prevents a supported T-cell program from being labelled NK. Unsupported or ambiguous assignments remain Uncertain. Plasma/B candidates are summarized as B_cell; Epithelial does not establish malignancy.

Then apply the fixed cell-level rare gates at expression ≥0.5:

- Mast: at least two of TPSAB1/TPSB2/CPA3 plus KIT or MS4A2.
- pDC: at least three of CLEC4C/LILRA4/IL3RA/TCF4/GZMB, including CLEC4C or LILRA4 and IL3RA or TCF4.
- Conflicting rare programs, or rare evidence together with at least two CD3 genes or both EPCAM/KRT19, yield Uncertain.

Cells without sufficient rare-marker evidence retain their cluster label. Thresholds 0.25 and 1.0 are sensitivity checks, not a search for maximum reference agreement. Marker rules are unchanged from the preceding development analysis and informed by [breast atlas stromal distinctions](https://www.nature.com/articles/s41588-021-00911-1), [pDC profiling](https://pmc.ncbi.nlm.nih.gov/articles/PMC10470411/) and [immune signatures](https://doi.org/10.3389/fimmu.2023.1194745). Reference labels, patient IDs and treatment metadata do not enter these annotation functions. This larger sample overlaps the source cohort used during rule development and is not an independent validation cohort.

| threshold | lineage | cells | jaccard_vs_primary |
| --- | --- | --- | --- |
| 0.25 | Mast_cell | 6 | 1.0 |
| 0.25 | pDC | 55 | 1.0 |
| 0.5 | Mast_cell | 6 | 1.0 |
| 0.5 | pDC | 55 | 1.0 |
| 1.0 | Mast_cell | 6 | 1.0 |
| 1.0 | pDC | 56 | 0.982 |


## Inspect stromal marker evidence

```{image} ../assets/tutorial-10000/rare_refinement-stromal_marker_review.png
:class: compact-scientific-plot
:align: center
:alt: Stromal-marker expression and detection across the newly fitted clusters. Use this evidence to review the unresolved endothelial/perivascular assignment.
```

Stromal-marker expression and detection across the newly fitted clusters. Use this evidence to review the unresolved endothelial/perivascular assignment.



## Outputs to use next

The optional diagnostic object is `rare_refinement/annotated_10000_refined.h5ad` under your output directory, with labels in `scknife_breast_refined`. Keep `fitted_subset.h5ad` and `activity_scores.parquet` alongside it. The main analysis does not load this diagnostic object; it continues to use the original `cellType` labels in `fitted_subset.h5ad`.

## View the generated plots

The code saves PNGs using a non-interactive plotting backend. Open the files in your output directory, or display one in Jupyter after the analysis completes:

```python
import os
from pathlib import Path
from IPython.display import Image, display
result_dir = Path(os.environ.get("SCKNIFE_TUTORIAL_OUTPUT", "tutorial/outputs_10000"))
display(Image(filename=str(result_dir / "rare_refinement/rare_refinement_umap.png")))
```

[Tutorial overview](../tutorial-10000.md) · [Previous: Set up and fit](10000-setup.md) · [Return to main tutorial](../tutorial-10000.md)
