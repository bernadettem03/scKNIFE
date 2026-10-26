import numpy as np
from sklearn.decomposition import NMF


def run_nmf(X, Rg, k=150):
    overlapping_idx = np.where(Rg.getnnz(axis=1) == 1)[0]
    X_overlap = X[overlapping_idx, :]
    nmf = NMF(n_components=k, init='nndsvda', max_iter=1000, random_state=42)
    W_overlap = nmf.fit_transform(X_overlap)
    H = nmf.components_
    return W_overlap, H, overlapping_idx


def initialize_W(W_overlap, Rg, overlapping_idx, n_nodes, k):
    W = np.zeros((n_nodes, k))
    node_idx = Rg[overlapping_idx].nonzero()[1]
    for i, node in enumerate(node_idx):
        W[node, :] = W_overlap[i, :]
    return W, node_idx
