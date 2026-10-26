import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import cg


def propagate_embeddings(W, L, labeled_idx, n_nodes, k, epsilon=1e-5):
    all_idx = np.arange(n_nodes)
    unlabeled_idx = all_idx[~np.isin(all_idx, labeled_idx)]
    L = L.tocsr()
    L_UU = L[unlabeled_idx, :][:, unlabeled_idx]
    L_UF = L[unlabeled_idx, :][:, labeled_idx]

    if not len(unlabeled_idx):
        return W, unlabeled_idx

    identity = epsilon * diags([1.0], [0], shape=L_UU.shape)
    b_all = -L_UF @ W[labeled_idx, :]
    W_U = np.zeros_like(b_all)

    for j in range(W.shape[1]):
        b = b_all[:, j]
        W_U[:, j], info = cg(L_UU + identity, b, rtol=1e-4, maxiter=1000)
        if info != 0:
            raise RuntimeError(f"Graph propagation did not converge for factor {j} (info={info})")
        W_U[:, j] = np.maximum(W_U[:, j], 0)
    for i, idx in enumerate(unlabeled_idx):
        W[idx, :] = W_U[i, :]
    return W, unlabeled_idx
