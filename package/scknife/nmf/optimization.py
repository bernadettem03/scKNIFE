import numpy as np
import anndata

from scknife.utils.helpers import store_H_to_obs


def run_graph_regularized_nmf(X, adata, W, Rg, A, D, L, k=150, lambda_reg=10,
                              beta_W=0.1, beta_H=0.1, n_iter=1000,
                              tol=1e-4, patience=10, random_state=42):
    print("Running scKNIFE Graph-Regularized NMF on CPU...")
    if n_iter < 1 or patience < 1:
        raise ValueError("n_iter and patience must be positive")
    if any(not np.isfinite(v) or v < 0 for v in (lambda_reg, beta_W, beta_H, tol)):
        raise ValueError("Penalties and tolerance must be finite and nonnegative")
    H = np.random.default_rng(random_state).random((k, X.shape[1]))
    epsilon = 1e-9
    loss_history = []
    prev_loss = None
    stale_epochs = 0
    X_norm_sq = np.linalg.norm(X, 'fro') ** 2

    for epoch in range(n_iter):
        RgW = Rg @ W
        numer_W = Rg.T @ (X @ H.T) + lambda_reg * A @ W
        denom_W = Rg.T @ (RgW @ H @ H.T) + lambda_reg * D @ W + beta_W + epsilon
        W *= numer_W / denom_W

        RgW = Rg @ W
        numer_H = RgW.T @ X
        denom_H = RgW.T @ RgW @ H + beta_H + epsilon
        H *= numer_H / denom_H

        recon_loss = np.linalg.norm(X - Rg @ W @ H, 'fro')**2
        smooth_loss = lambda_reg * np.trace(W.T @ L @ W)
        sparse_loss_W = beta_W * np.sum(W)
        sparse_loss_H = beta_H * np.sum(H)
        total_loss = 0.5 * recon_loss + 0.5 * smooth_loss + sparse_loss_W + sparse_loss_H
        total_loss_normalized = total_loss / X_norm_sq
        loss_history.append(total_loss_normalized)

        if prev_loss is not None:
            rel_change = abs(prev_loss - total_loss_normalized) / max(abs(prev_loss), 1e-10)
            if rel_change < tol:
                stale_epochs += 1
                if stale_epochs >= patience:
                    print(f"Epoch {epoch}: converged (rel_change={rel_change:.6f} < {tol})")
                    break
            else:
                stale_epochs = 0

        if epoch % 100 == 0 or epoch == n_iter - 1:
            print(f"Epoch {epoch}: Normalized Loss = {total_loss_normalized:.6f}")

        prev_loss = total_loss_normalized

    adata.uns['W_graph'] = W
    store_H_to_obs(adata, H, prefix="H_graph")

    print(f"Stored W_graph in `adata.uns['W_graph']` ({loss_history[-1]:.6f} final normalized loss)")
    print("Stored H_graph components in `adata.obs`")
    return W, H, loss_history
