# Hallmark-focused graph collection

This preserved GraphML collection provides a focused view of the biological network. Its name identifies the collection focus; the file retains relationships from the source labels listed below.

## File identity

- File: `hallmark.graphml`
- Nodes: 4,434
- Edges: 50,474
- SHA-256: `bacfb2a2fce777c0c4a7a1431db38bf94f250e4f3ba3af99ad2f921e8fc0aa0f`

## Relationship sources

Counts below are edge-source memberships. An edge with multiple source labels contributes once to each label, so these counts can sum to more than the number of edges.

| Recorded source label | Edge memberships | Resource reference |
|---|---:|---|
| `ARCHS4 TFs Coexp Edges` | 14,814 | [ARCHS4](https://doi.org/10.1038/s41467-018-03751-6) |
| `Human-GEM` | 26 | [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) |
| `breast-PPI` | 28,206 | Source identifier retained from the supplied graph |
| `cellchat` | 669 | [CellChat](https://doi.org/10.1038/s41467-021-21246-9) |
| `hallmark` | 7,322 | [Hallmark / MSigDB](https://doi.org/10.1016/j.cels.2015.12.004) |
| `rxn_gene` | 11 | [Metabolic task/reaction resources](https://github.com/earmingol/scCellFie) |

## Use and citation

Load the complete collection with `networkx.read_graphml()`. Select this graph before tuning and keep its node order and file identity with the analysis. For source-specific analyses, use the stored edge `source` attribute to identify the relationships of interest.

Cite scKNIFE and the biological resources used in the analysis. See the [resource references](../../../docs/graph-resources.md) and [graph selection guide](../../../docs/graphs.md). Resource-specific terms accompany the biological data; the scKNIFE software license applies to the software.
