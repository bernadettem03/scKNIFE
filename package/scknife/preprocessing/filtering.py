import numpy as np
import pandas as pd


def compute_mean_expression(adata):
    X = adata.X
    if hasattr(X, "toarray"):
        mean_expression = np.asarray(X.mean(axis=0)).ravel()
    else:
        mean_expression = X.mean(axis=0)
    return pd.Series(mean_expression, index=adata.var_names, name='E_score')


def filter_top_genes_by_mean_expression(adata, top_n=10000):
    E_scores = compute_mean_expression(adata)
    top_genes = E_scores.sort_values(ascending=False).head(top_n).index
    return adata[:, top_genes].copy()
