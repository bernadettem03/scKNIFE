import anndata as ad,numpy as np,pandas as pd,pytest
from scknife.celltype.rare_refinement import refine_rare_breast_labels

def test_combinations_not_single_markers_and_conflicts():
    genes=['TPSAB1','TPSB2','CPA3','KIT','CLEC4C','IL3RA','TCF4','GZMB','CD3D','CD3E']
    x=np.zeros((5,len(genes)))
    x[0,[0,1,3]]=2
    x[1,[4,5,6,7]]=2
    x[2,7]=3
    x[3,[4,5,6,8,9]]=2
    x[4,[0,1,3,4,5,6]]=2
    a=ad.AnnData(x,var=pd.DataFrame(index=genes));base=pd.Series('Myeloid',index=a.obs_names)
    result,evidence,_=refine_rare_breast_labels(a,base)
    assert result.tolist()==['Mast_cell','pDC','Myeloid','Uncertain','Uncertain']
    a.obs['cellType']='deliberately wrong'
    again,*_=refine_rare_breast_labels(a,base.iloc[::-1])
    pd.testing.assert_series_equal(result,again)

def test_missing_anchor_does_not_rescue_and_bad_input():
    a=ad.AnnData(np.ones((2,3))*3,var=pd.DataFrame(index=['TCF4','GZMB','IL3RA']))
    base=pd.Series('Myeloid',index=a.obs_names)
    result,*_=refine_rare_breast_labels(a,base)
    assert result.eq('Myeloid').all()
    with pytest.raises(ValueError):refine_rare_breast_labels(a,base,min_expression=0)
