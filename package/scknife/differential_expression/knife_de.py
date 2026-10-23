import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests


def cohens_d(x, y):
    if len(x) < 2 or len(y) < 2:
        return np.nan
    s1 = np.var(x, ddof=1)
    s2 = np.var(y, ddof=1)
    sp = np.sqrt((s1 + s2) / 2)
    if sp == 0:
        return 0
    return (np.mean(y) - np.mean(x)) / sp


def run_scanpy_KNIFE_DE(
    mt_df,
    meta_df,
    cell_type_key,
    condition_key,
    contrasts,
    features=None,
    min_cells=0,
):
    """
    Scanpy Wilcoxon DE with KNIFE logic:
    cell-type aware, effect-size aware, globally FDR corrected.
    """
    adata = ad.AnnData(
        X=mt_df.values,
        obs=meta_df.copy(),
        var=pd.DataFrame(index=mt_df.columns),
    )

    adata.obs[cell_type_key] = adata.obs[cell_type_key].astype("category")
    adata.obs[condition_key] = adata.obs[condition_key].astype("category")

    if features is None:
        features = mt_df.columns.tolist()

    all_results = []

    for ct in adata.obs[cell_type_key].cat.categories:
        ad_ct = adata[adata.obs[cell_type_key] == ct].copy()

        for c1, c2 in contrasts:
            mask1 = ad_ct.obs[condition_key] == c1
            mask2 = ad_ct.obs[condition_key] == c2

            if mask1.sum() == 0 or mask2.sum() == 0 or min(mask1.sum(), mask2.sum()) < min_cells:
                continue

            ad_pair = ad_ct[mask1 | mask2].copy()

            sc.tl.rank_genes_groups(
                ad_pair,
                groupby=condition_key,
                groups=[c2],
                reference=c1,
                method="wilcoxon",
            )

            res = sc.get.rank_genes_groups_df(ad_pair, group=c2)

            res = res[res["names"].isin(features)]

            for _, row in res.iterrows():
                gene = row["names"]

                x = ad_pair[ad_pair.obs[condition_key] == c1, gene].X
                y = ad_pair[ad_pair.obs[condition_key] == c2, gene].X

                if issparse(x): x = x.A.flatten()
                else: x = np.asarray(x).flatten()

                if issparse(y): y = y.A.flatten()
                else: y = np.asarray(y).flatten()

                m1, m2 = np.mean(x), np.mean(y)
                log2fc = np.log2(m2 / m1) if m1 > 0 else np.nan

                all_results.append({
                    "cell_type": ct,
                    "feature": gene,
                    "group1": c1,
                    "group2": c2,
                    "log2FC": log2fc,
                    "test_statistic": row["scores"],
                    "p_value": row["pvals"],
                    "cohens_d": cohens_d(x, y),
                    "n_group1": len(x),
                    "n_group2": len(y),
                    "median_group1": np.median(x),
                    "median_group2": np.median(y),
                    "median_diff": np.median(y) - np.median(x),
                })

    df = pd.DataFrame(all_results)

    if not df.empty:
        df["adj_p_value"] = multipletests(df["p_value"], method="fdr_bh")[1]

    return df
