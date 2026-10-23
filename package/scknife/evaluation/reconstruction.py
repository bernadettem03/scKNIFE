import numpy as np
from scipy.stats import pearsonr
from sklearn.decomposition import NMF

from scknife.utils.helpers import store_H_to_obs


def evaluate_reconstruction(X, Rg, W, H):
    X_recon = Rg @ W @ H
    corr_flat, _ = pearsonr(X.flatten(), X_recon.flatten())
    fro_error = np.linalg.norm(X - X_recon, ord='fro')
    rel_error = fro_error / np.linalg.norm(X, ord='fro')
    return corr_flat, fro_error, rel_error


def run_standard_nmf_and_evaluate(X, adata, k=150):
    print("Running Standard NMF...")
    nmf = NMF(n_components=k, init='random', max_iter=1000, random_state=42)
    W_std = nmf.fit_transform(X)
    H_std = nmf.components_
    X_hat_std = W_std @ H_std

    corr = np.corrcoef(X.flatten(), X_hat_std.flatten())[0, 1]
    fro_error = np.linalg.norm(X - X_hat_std, ord='fro')
    rel_error = fro_error / np.linalg.norm(X, ord='fro')

    adata.uns['W_std'] = W_std
    store_H_to_obs(adata, H_std, prefix="H_std")

    print("Stored W_std in `adata.uns['W_std']`")
    print("Stored H_std components in `adata.obs`")
    print(f"Std NMF -> Corr: {corr:.4f}, Frobenius: {fro_error:.2f}, Relative: {rel_error:.4f}")

    return W_std, H_std, X_hat_std, {
        "correlation": corr,
        "frobenius_error": fro_error,
        "relative_error": rel_error
    }


def compare_nmf_models(X, Rg, W_graph, H_graph, W_std, H_std, adata=None):
    X_scknife_recon = Rg @ W_graph @ H_graph
    X_std_recon = W_std @ H_std

    scknife_corr = np.corrcoef(X.flatten(), X_scknife_recon.flatten())[0, 1]
    scknife_fro = np.linalg.norm(X - X_scknife_recon, 'fro')
    scknife_rel = scknife_fro / np.linalg.norm(X, 'fro')

    std_nmf_corr = np.corrcoef(X.flatten(), X_std_recon.flatten())[0, 1]
    std_nmf_fro = np.linalg.norm(X - X_std_recon, 'fro')
    std_nmf_rel = std_nmf_fro / np.linalg.norm(X, 'fro')

    metrics = {
        "scknife_corr": float(scknife_corr),
        "scknife_fro": float(scknife_fro),
        "scknife_rel": float(scknife_rel),
        "std_nmf_corr": float(std_nmf_corr),
        "std_nmf_fro": float(std_nmf_fro),
        "std_nmf_rel": float(std_nmf_rel),
    }

    if adata is not None:
        adata.uns['recon_metrics'] = metrics

    print("\nStandard NMF:")
    print(f"   Pearson Correlation     : {std_nmf_corr:.4f}")
    print(f"   Frobenius Error         : {std_nmf_fro:.4f}")
    print(f"   Relative Frobenius Error: {std_nmf_rel:.4f}")

    print("\nscKNIFE (Graph-Regularized NMF):")
    print(f"   Pearson Correlation     : {scknife_corr:.4f}")
    print(f"   Frobenius Error         : {scknife_fro:.4f}")
    print(f"   Relative Frobenius Error: {scknife_rel:.4f}")
