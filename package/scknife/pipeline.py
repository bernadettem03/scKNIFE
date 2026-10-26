from scknife.preprocessing.gene_mapping import extract_expression_matrix
from scknife.graph.laplacian import build_graph_laplacian
from scknife.graph.mappings import build_Rg_mapping
from scknife.nmf.initialization import run_nmf, initialize_W
from scknife.nmf.propagation import propagate_embeddings
from scknife.nmf.optimization import run_graph_regularized_nmf
from scknife.plotting.losses import plot_losses
from scknife.evaluation.reconstruction import run_standard_nmf_and_evaluate, compare_nmf_models


def run_full_pipeline(
    adata,
    G,
    k=30,
    lambda_reg=10,
    beta_W=0.1,
    beta_H=0.1,
    n_iter=500,
    verbose=True,
    tol=1e-4,
    random_state=42,
    max_dense_bytes=2_000_000_000,
):
    print("Starting Full scKNIFE (Graph-Regularized NMF) Pipeline...")
    print(f"Params - k={k}, lambda={lambda_reg}, beta_W={beta_W}, beta_H={beta_H}, Iterations={n_iter}\n")

    X, adata = extract_expression_matrix(adata, G, max_dense_bytes=max_dense_bytes)
    n_genes, n_cells = X.shape
    if not isinstance(k, int) or not 1 <= k <= min(X.shape):
        raise ValueError("k must be between 1 and min(n_shared_genes, n_cells)")
    if verbose:
        print(f"Genes retained: {n_genes}")
        print(f"Cells retained: {n_cells}\n")

    A, D, L = build_graph_laplacian(G)
    if verbose:
        print(f"Graph nodes: {A.shape[0]}")
        print(f"Graph edges: {A.nnz // 2}\n")

    Rg, mapping_idx, node_list = build_Rg_mapping(adata, G)
    if verbose:
        print(f"Rg shape: {Rg.shape}")
        print(f"Genes mapped to KG nodes: {len(mapping_idx)}\n")

    W_overlap, H_overlap, overlapping_idx = run_nmf(X, Rg, k=k)
    if verbose:
        print(f"NMF on {len(overlapping_idx)} overlapping genes -> W_overlap: {W_overlap.shape}, H: {H_overlap.shape}\n")

    W_init, labeled_idx = initialize_W(W_overlap, Rg, overlapping_idx, n_nodes=Rg.shape[1], k=k)
    if verbose:
        print(f"Initialized W: {W_init.shape}, Labeled nodes: {len(labeled_idx)}\n")

    W_propagated, unlabeled_idx = propagate_embeddings(W_init, L, labeled_idx, n_nodes=Rg.shape[1], k=k)
    if verbose:
        print(f"Propagated W to unlabeled nodes: {len(unlabeled_idx)}\n")

    W_graph, H_graph, loss_history = run_graph_regularized_nmf(
        X, adata, W_propagated, Rg, A, D, L, k=k,
        lambda_reg=lambda_reg, beta_W=beta_W, beta_H=beta_H, n_iter=n_iter,
        tol=tol, random_state=random_state,
    )

    # Store loss history for downstream evaluation (joint grid, tuning)
    adata.uns['loss_history'] = loss_history
    plot_losses(loss_history)

    W_std, H_std, X_hat_std, std_metrics = run_standard_nmf_and_evaluate(X, adata, k=k)

    compare_nmf_models(X, Rg, W_graph, H_graph, W_std, H_std, adata=adata)
    adata.uns['Rg'] = Rg

    print("\nPipeline complete!")
    print("Results stored in `adata.obs` and `adata.uns`:\n"
          "   - `W_std`, `W_graph`\n"
          "   - `H_std_{i}`, `H_graph_{i}`\n"
          "   - `loss_history`\n"
          "   - `Rg`, `node_names`\n")

    return adata
