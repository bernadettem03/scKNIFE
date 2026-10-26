# Reaction–gene graph collection

This preserved GraphML collection provides a focused view of the biological network. Its name identifies the collection focus; the file retains relationships from the source labels listed below.

## File identity

- File: `rxn_by_gene.graphml`
- Nodes: 1,705
- Edges: 5,992
- SHA-256: `3793d707f6adca68eccd5e65f20b386c4be7cfdbd03a6a61a4d5a98d12a85866`

## Relationship sources

Counts below are edge-source memberships. An edge with multiple source labels contributes once to each label, so these counts can sum to more than the number of edges.

| Recorded source label | Edge memberships | Resource reference |
|---|---:|---|
| `ARCHS4 TFs Coexp Edges` | 1,141 | [ARCHS4](https://doi.org/10.1038/s41467-018-03751-6) |
| `Human-GEM` | 124 | [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) |
| `breast-PPI` | 2,679 | Source identifier retained from the supplied graph |
| `rxn_gene` | 2,315 | [Metabolic task/reaction resources](https://github.com/earmingol/scCellFie) |

## Use and citation

Load the complete collection with `networkx.read_graphml()`. Select this graph before tuning and keep its node order and file identity with the analysis. For source-specific analyses, use the stored edge `source` attribute to identify the relationships of interest.

Cite scKNIFE and the biological resources used in the analysis. See the [resource references](../../../docs/graph-resources.md) and [graph selection guide](../../../docs/graphs.md). Resource-specific terms accompany the biological data; the scKNIFE software license applies to the software.
