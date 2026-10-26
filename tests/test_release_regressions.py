import numpy as np
import pandas as pd
import pytest
import anndata as ad
import networkx as nx
from scipy import sparse
from scknife import compute_activity_scores, run_full_pipeline, run_cellstate_assignment
from scknife.utils.helpers import store_H_to_obs
from scknife.preprocessing.gene_mapping import extract_expression_matrix
from scknife.signaling import activity_lr as lr


def test_activity_numeric_order_and_invalid_factors():
    a = ad.AnnData(np.ones((1, 1)))
    a.obs['H_graph_1'] = [20.]
    a.obs['H_graph_0'] = [10.]
    a.uns.update(W_graph=np.array([[1., 2.]]), node_names=['node'])
    compute_activity_scores(a)
    assert a.uns['activity_scores'].iloc[0, 0] == 50
    a.obs = a.obs.drop(columns='H_graph_0')
    with pytest.raises(ValueError, match='exactly once'):
        compute_activity_scores(a)


def test_refit_replaces_old_factors():
    a = ad.AnnData(np.ones((3, 2)))
    store_H_to_obs(a, np.ones((4, 3)), 'H_graph')
    store_H_to_obs(a, np.ones((2, 3))*2, 'H_graph')
    assert list(a.obs) == ['H_graph_0', 'H_graph_1']
    assert (a.obs.to_numpy() == 2).all()


def test_sparse_dense_seeded_pipeline_and_guard():
    a = ad.AnnData(np.random.default_rng(3).random((12, 5)),
                   var=pd.DataFrame(index=list('ABCDE')))
    g = nx.path_graph(list('ABCDEF'))
    b = a.copy(); b.X = sparse.csr_matrix(b.X)
    with pytest.raises(ValueError, match='Dense expression'):
        extract_expression_matrix(b, g, max_dense_bytes=1)
    args = dict(k=2, n_iter=4, tol=0, random_state=7, verbose=False)
    x = run_full_pipeline(a, g, **args)
    y = run_full_pipeline(b, g, **args)
    np.testing.assert_allclose(x.uns['W_graph'], y.uns['W_graph'])
    np.testing.assert_allclose(x.obs[['H_graph_0','H_graph_1']], y.obs[['H_graph_0','H_graph_1']])
    assert np.all(np.diff(x.uns['loss_history']) <= 1e-8)
    with pytest.raises(ValueError, match='nonnegative'):
        bad = a.copy(); bad.X[0,0] = -1
        extract_expression_matrix(bad, g)


def test_cellstate_custom_cluster_key():
    a = ad.AnnData(np.ones((4, 1)))
    a.obs['custom'] = ['alpha', 'alpha', 'beta', 'beta']
    activity = pd.DataFrame([[3,3,0,0],[0,0,4,4]], index=['s1','s2'], columns=a.obs_names)
    result, _, _ = run_cellstate_assignment(a, activity, ['s1','s2'], cluster_key='custom')
    assert result.obs.bgnmf_cell_state.tolist() == ['s1','s1','s2','s2']


def test_permutation_resume_and_input_invalidation(tmp_path, monkeypatch):
    monkeypatch.setattr(lr, 'DATA_DIR', tmp_path)
    activity = np.random.default_rng(4).random((8,2))
    labels = np.repeat([0,1],4)
    pairs = pd.DataFrame(dict(ligand_index=[0],receptor_index=[1],confidence=[1.],cofactor_weight=[1.]))
    specificity = np.ones((2,2))
    observed = np.full((2,2), 1e9)  # no null exceedances
    def run(n, force=False):
        return lr.run_permutations(activity, labels, observed, specificity, pairs, 2,n,2,force,'same-labels')[0]
    np.testing.assert_allclose(run(3), 1/4)
    resumed = run(7)
    np.testing.assert_allclose(resumed, run(7, True))
    with np.load(tmp_path/'permutation_checkpoint.npz') as saved:
        before = str(saved['input_signature'])
    activity[0,0] += 5
    changed = run(7)
    with np.load(tmp_path/'permutation_checkpoint.npz') as saved:
        assert before != str(saved['input_signature'])
    np.testing.assert_allclose(changed, run(7, True))
    assert (changed > 0).all()
