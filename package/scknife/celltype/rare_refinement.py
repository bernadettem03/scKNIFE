"""Conservative cell-level rare-lineage refinement of existing breast labels."""
import numpy as np
import pandas as pd
from scipy import sparse

RARE_PANELS = {
    'Mast_cell': ('TPSAB1','TPSB2','CPA3','KIT','MS4A2'),
    'pDC': ('CLEC4C','LILRA4','IL3RA','TCF4','GZMB'),
}

def refine_rare_breast_labels(expression, base_labels, *, min_expression=0.5):
    """Refine with coexpressed markers on normalize_total(1e4)/log1p expression.

    Mast: >=2 of TPSAB1/TPSB2/CPA3 plus KIT or MS4A2.
    pDC: >=3 panel genes, including CLEC4C or LILRA4 AND IL3RA or TCF4.
    Every gate uses expression >= min_expression. Two CD3 markers or both
    EPCAM/KRT19 conflict with a rare call; simultaneous rare programs also
    yield Uncertain. Otherwise preserve the input label. No reference metadata
    or cell identities influence scoring. This heuristic is not a probability.
    """
    if min_expression<=0:raise ValueError('min_expression must be positive')
    if not expression.obs_names.is_unique or not expression.var_names.is_unique:
        raise ValueError('Expression identifiers must be unique')
    if not isinstance(base_labels,pd.Series) or not base_labels.index.is_unique:
        raise ValueError('Base labels must be a uniquely indexed Series')
    labels=base_labels.reindex(expression.obs_names)
    if labels.isna().any():raise ValueError('Base labels must cover all cells')
    requested=sorted(set(sum(RARE_PANELS.values(),()))|{'CD3D','CD3E','CD3G','EPCAM','KRT19'})
    present=[g for g in requested if g in expression.var_names]
    x=expression[:,present].X;x=x.toarray() if sparse.issparse(x) else np.asarray(x)
    if not np.isfinite(x).all() or (x<0).any():raise ValueError('Expression must be finite and nonnegative')
    measured=pd.DataFrame(x,index=expression.obs_names,columns=present)
    hits=measured.ge(min_expression).reindex(columns=requested,fill_value=False)
    mast=(hits[['TPSAB1','TPSB2','CPA3']].sum(axis=1)>=2)&hits[['KIT','MS4A2']].any(axis=1)
    pdc=(hits[list(RARE_PANELS['pDC'])].sum(axis=1)>=3)&hits[['CLEC4C','LILRA4']].any(axis=1)&hits[['IL3RA','TCF4']].any(axis=1)
    conflict=(hits[['CD3D','CD3E','CD3G']].sum(axis=1)>=2)|hits[['EPCAM','KRT19']].all(axis=1)
    ambiguous=(mast&pdc)|((mast|pdc)&conflict)
    refined=labels.astype(str).copy();refined.name='scknife_breast_refined'
    refined.loc[mast&~ambiguous]='Mast_cell';refined.loc[pdc&~ambiguous]='pDC';refined.loc[ambiguous]='Uncertain'
    evidence=pd.DataFrame({'base_label':labels,'refined_label':refined,'mast_evidence':mast,'pdc_evidence':pdc,'conflicting_lineage':conflict,'ambiguous_rare_evidence':ambiguous},index=expression.obs_names)
    evidence['changed']=evidence.base_label.ne(evidence.refined_label)
    coverage=pd.DataFrame({'gene':requested,'available':[g in present for g in requested]})
    return refined,evidence,coverage
