# CellChat-focused graph collection

This preserved GraphML collection provides a focused view of the biological network. Its name identifies the collection focus; the file retains relationships from the source labels listed below.

## File identity

- File: `cellchat.graphml`
- Nodes: 1,160
- Edges: 4,144
- SHA-256: `3aa83c8a59914256bac310e414a41ec2b35d0907df899fbe9aeda0b7c86b7d8a`

## Relationship sources

Counts below are edge-source memberships. An edge with multiple source labels contributes once to each label, so these counts can sum to more than the number of edges.

| Recorded source label | Edge memberships | Resource reference |
|---|---:|---|
| `ARCHS4 TFs Coexp Edges` | 1,349 | [ARCHS4](https://doi.org/10.1038/s41467-018-03751-6) |
| `breast-PPI` | 633 | Source identifier retained from the supplied graph |
| `cellchat` | 2,287 | [CellChat](https://doi.org/10.1038/s41467-021-21246-9) |

## Use and citation

Load the complete collection with `networkx.read_graphml()`. Select this graph before tuning and keep its node order and file identity with the analysis. For source-specific analyses, use the stored edge `source` attribute to identify the relationships of interest.

Cite scKNIFE and the biological resources used in the analysis. See the [resource references](../../../docs/graph-resources.md) and [graph selection guide](../../../docs/graphs.md). Resource-specific terms accompany the biological data; the scKNIFE software license applies to the software.
