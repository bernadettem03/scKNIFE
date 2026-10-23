# Marker-focused graph collection

This preserved GraphML collection provides a focused view of the biological network. Its name identifies the collection focus; the file retains relationships from the source labels listed below.

## File identity

- File: `markers.graphml`
- Nodes: 1,064
- Edges: 3,745
- SHA-256: `6f1a9f4bc19ff5e91c221d975e1628fadda392f5674370f4de44220cd0b7ec71`

## Relationship sources

Counts below are edge-source memberships. An edge with multiple source labels contributes once to each label, so these counts can sum to more than the number of edges.

| Recorded source label | Edge memberships | Resource reference |
|---|---:|---|
| `ARCHS4 TFs Coexp Edges` | 1,867 | [ARCHS4](https://doi.org/10.1038/s41467-018-03751-6) |
| `breast-PPI` | 591 | Source identifier retained from the supplied graph |
| `cellchat` | 90 | [CellChat](https://doi.org/10.1038/s41467-021-21246-9) |
| `marker_graph` | 1,230 | [Marker resources](../../../docs/graph-resources.md) |

## Use and citation

Load the complete collection with `networkx.read_graphml()`. Select this graph before tuning and keep its node order and file identity with the analysis. For source-specific analyses, use the stored edge `source` attribute to identify the relationships of interest.

Cite scKNIFE and the biological resources used in the analysis. See the [resource references](../../../docs/graph-resources.md) and [graph selection guide](../../../docs/graphs.md). Resource-specific terms accompany the biological data; the scKNIFE software license applies to the software.
