"""Breast-context broad-lineage annotation from measured marker expression.

This is an explicit expression-based alternative to graph-node argmax labels.
Reference annotations and clinical metadata are not inputs to the scoring rule.
"""
import numpy as np
import pandas as pd
from scipy import sparse

BREAST_MARKERS = {
    'Epithelial': ('EPCAM','KRT8','KRT18','KRT19','CDH1'),
    'Fibroblast': ('DCN','LUM','COL1A1','COL1A2','COL3A1','PDGFRA'),
    'Endothelial': ('PECAM1','VWF','CLDN5','EMCN','KDR'),
    'Perivascular': ('RGS5','CSPG4','MCAM','NOTCH3','MYH11'),
    'T_cell': ('CD3D','CD3E','CD3G','TRAC'),
    'NK_cell': ('GNLY','KLRD1','NCR1','PRF1'),
    'B_cell': ('MS4A1','CD79A','CD79B','CD19','CD22'),
    'Plasma_cell': ('MZB1','JCHAIN','SDC1','DERL3'),
    'Myeloid': ('LYZ','LST1','TYROBP','FCER1G','AIF1','CTSS'),
    'pDC': ('CLEC4C','IL3RA','TCF4','GZMB'),
    'Mast_cell': ('TPSAB1','TPSB2','CPA3','KIT','MS4A2'),
}

def run_breast_celltype_assignment(expression, cluster_labels, *, min_score=0.5,
        min_margin=0.15, min_detection=0.1, min_markers=2):
    """Annotate clusters using normalized/log1p measured expression.

    Standardize each marker's cluster mean ACROSS clusters (population SD),
    clip z-scores to [-3, 3], then average available markers per lineage.
    Require >=min_markers detected in >=min_detection fraction of the cluster,
    a positive score >=min_score, and a top-vs-next-eligible-lineage score gap >=min_margin.
    Unsupported panels are excluded before ranking. CD3 support excludes the
    NK candidate; plasma and B candidates share the broad B_cell output label.
    Scores/margins are heuristic evidence, not calibrated probabilities.
    Missing/weak/ambiguous evidence yields Uncertain. Returns cell labels,
    cluster decisions, lineage scores, marker detection, and marker coverage.
    """
    if not expression.obs_names.is_unique or not expression.var_names.is_unique:
        raise ValueError('Cell and gene identifiers must be unique')
    if not isinstance(cluster_labels,pd.Series) or not cluster_labels.index.is_unique:
        raise ValueError('Supply uniquely indexed cluster labels')
    if not expression.obs_names.isin(cluster_labels.index).all():
        raise ValueError('Cluster labels must cover every expression cell')
    labels=cluster_labels.reindex(expression.obs_names)
    if labels.isna().any() or labels.nunique()<2:
        raise ValueError('Need at least two non-missing clusters')
    if min_markers<1 or not 0<=min_detection<=1 or min_margin<0 or min_score<0:
        raise ValueError('Invalid annotation thresholds')
    genes=sorted(set(g for panel in BREAST_MARKERS.values() for g in panel)&set(expression.var_names))
    if not genes:raise ValueError('No breast marker genes present')
    x=expression[:,genes].X
    x=x.toarray() if sparse.issparse(x) else np.asarray(x)
    if not np.isfinite(x).all() or (x<0).any():raise ValueError('Expression must be finite and nonnegative')
    frame=pd.DataFrame(x,index=expression.obs_names,columns=genes)
    means=frame.groupby(labels,observed=True,sort=True).mean()
    detection=frame.gt(0).groupby(labels,observed=True,sort=True).mean()
    std=means.std(axis=0,ddof=0)
    z=means.sub(means.mean(axis=0)).div(std.where(std>1e-12,1)).clip(-3,3)
    scores=pd.DataFrame(index=means.index); support=pd.DataFrame(index=means.index); coverage=[]
    for lineage,panel in BREAST_MARKERS.items():
        present=[g for g in panel if g in genes]
        coverage.append({'lineage':lineage,'available':','.join(present),'missing':','.join(g for g in panel if g not in genes)})
        scores[lineage]=z[present].mean(axis=1) if len(present)>=min_markers else -np.inf
        support[lineage]=(detection[present]>=min_detection).sum(axis=1)
    # Exclude unsupported candidates BEFORE ranking. A rare marker's high z-score
    # cannot win when too few cells detect the panel.
    eligible=scores.where(support>=min_markers,-np.inf)
    # T-cell receptor/CD3 evidence distinguishes cytotoxic T cells from NK cells.
    cd3=[g for g in ('CD3D','CD3E','CD3G') if g in detection]
    t_supported=(detection[cd3]>=min_detection).sum(axis=1)>=2
    eligible.loc[t_supported,'NK_cell']=-np.inf
    # Report B/plasma at broad B-lineage granularity; preserve candidate subtype.
    family={n:('B_cell' if n=='Plasma_cell' else n) for n in BREAST_MARKERS}
    rows=[]
    for cluster,row in eligible.iterrows():
        ranked=row.sort_values(ascending=False,kind='stable')
        top=ranked.index[0];best=float(ranked.iloc[0])
        alternatives=[float(v) for n,v in ranked.items() if family[n]!=family[top] and np.isfinite(v)]
        second=max(alternatives) if alternatives else 0.0
        margin=best-second if np.isfinite(best) else 0.0
        reason='accepted'
        if not np.isfinite(best):reason='insufficient marker support'
        elif best<min_score:reason='weak enrichment'
        elif margin<min_margin:reason='ambiguous lineage scores'
        rows.append({'cluster':cluster,'candidate':top,'label':family[top] if reason=='accepted' else 'Uncertain',
                     'score':best,'margin':margin,'detected_markers':int(support.loc[cluster,top]),'reason':reason,
                     't_cell_receptor_supported':bool(t_supported.loc[cluster]),
                     'cells':int(labels.eq(cluster).sum())})
    decisions=pd.DataFrame(rows).set_index('cluster')
    cell_labels=labels.map(decisions.label).rename('scknife_breast_cell_type')
    return cell_labels,decisions,scores,detection,pd.DataFrame(coverage)
