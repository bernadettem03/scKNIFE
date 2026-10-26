import numpy as np,pandas as pd,anndata as ad,pytest
from scknife.celltype.breast import run_breast_celltype_assignment

def test_marker_lineages_and_reference_independence():
    genes=['DCN','LUM','CD3D','CD3E']; x=np.array([[5,4,0,0],[4,5,0,0],[0,0,5,4],[0,0,4,5]])
    a=ad.AnnData(x,obs=pd.DataFrame(index=['a','b','c','d']),var=pd.DataFrame(index=genes))
    clusters=pd.Series([0,0,1,1],index=a.obs_names)
    labels,*_=run_breast_celltype_assignment(a,clusters)
    assert labels.tolist()==['Fibroblast','Fibroblast','T_cell','T_cell']
    a.obs['cellType']=['wrong']*4
    shuffled,*_=run_breast_celltype_assignment(a,clusters.iloc[::-1])
    pd.testing.assert_series_equal(labels,shuffled)

def test_weak_or_ambiguous_evidence_is_uncertain():
    a=ad.AnnData(np.ones((4,4)),var=pd.DataFrame(index=['DCN','LUM','CD3D','CD3E']))
    labels,*_=run_breast_celltype_assignment(a,pd.Series([0,0,1,1],index=a.obs_names))
    assert set(labels)=={'Uncertain'}
    a.X=np.array([[5,5,5,5],[5,5,5,5],[0,0,0,0],[0,0,0,0]])
    labels,decisions,*_=run_breast_celltype_assignment(a,pd.Series([0,0,1,1],index=a.obs_names))
    assert set(labels)=={'Uncertain'}
    assert decisions.loc[0,'reason']=='ambiguous lineage scores'

def test_missing_markers_and_misalignment():
    a=ad.AnnData(np.ones((4,1)),var=pd.DataFrame(index=['DCN']))
    labels,*_=run_breast_celltype_assignment(a,pd.Series([0,0,1,1],index=a.obs_names))
    assert set(labels)=={'Uncertain'}
    with pytest.raises(ValueError,match='cover every'):
        run_breast_celltype_assignment(a,pd.Series([0,1],index=['x','y']))


def test_cd3_support_prevents_cytotoxic_t_as_nk():
    genes=['CD3D','CD3E','CD3G','TRAC','GNLY','KLRD1','NCR1','PRF1','DCN','LUM']
    a=ad.AnnData(np.array([[4,4,4,4,9,9,9,9,0,0],[4,4,4,4,9,9,9,9,0,0],[0,0,0,0,0,0,0,0,5,5],[0,0,0,0,0,0,0,0,5,5]]),var=pd.DataFrame(index=genes))
    labels,decisions,*_=run_breast_celltype_assignment(a,pd.Series([0,0,1,1],index=a.obs_names))
    assert decisions.loc[0,'t_cell_receptor_supported']
    assert labels.iloc[0]=='T_cell'
