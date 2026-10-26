---
orphan: true
---

# Concepts and data model

scKNIFE combines a cell-by-gene expression matrix with a biological graph containing genes and higher-level nodes such as reactions, metabolic tasks, pathways, marker programs, cell states, or ligand-receptor knowledge.

## Learned matrices

- **W** stores graph-node loadings across latent factors.
- **H** stores cell loadings across the same factors.
- **Rg** maps expression genes to graph nodes.
- Reconstructed graph-node activity can be summarized by cell type, state, patient, or condition.

## AnnData contract

| Location | Content |
|---|---|
| `adata.uns["W_graph"]` | graph-regularized node loadings |
| `adata.obs["H_graph_{i}"]` | per-cell graph-factor loadings |
| `adata.uns["loss_history"]` | optimization history |
| `adata.uns["Rg"]` | gene-to-node mapping |
| `adata.uns["node_names"]` | graph-node order |

Clustering adds H-space cluster labels to `adata.obs`; existing cell-type labels remain grouping metadata. UMAP and H-space representations are stored in `adata.obsm`.

## Model output versus metadata

Cell type, patient, site, and response labels do not generate scKNIFE factors. They are joined by cell identity after fitting and used to summarize model-derived activity. When authentic sample identifiers are available, metadata can group cell-level scores into sample summaries. The main 10,000-cell tutorial retains authentic patient IDs and timepoints and presents exploratory cell-level metabolic and GO tests alongside cell–cell signaling summaries.
