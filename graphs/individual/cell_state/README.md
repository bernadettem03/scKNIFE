# Cell-state graph collection

This preserved GraphML collection provides a focused view of the biological network. Its name identifies the collection focus; the file retains relationships from the source labels listed below.

## File identity

- File: `cell_state.graphml`
- Nodes: 498
- Edges: 696
- SHA-256: `0577fe985d57f290a8bcad8199da73a63a336b9270b7a66c59d879975deab3bc`

## Relationship sources

Counts below are edge-source memberships. An edge with multiple source labels contributes once to each label, so these counts can sum to more than the number of edges.

| Recorded source label | Edge memberships | Resource reference |
|---|---:|---|
| `cellchat` | 43 | [CellChat](https://doi.org/10.1038/s41467-021-21246-9) |
| `cellmarker2.0` | 653 | [CellMarker 2.0](https://doi.org/10.1093/nar/gkac947) |

## Use and citation

Load the complete collection with `networkx.read_graphml()`. Select this graph before tuning and keep its node order and file identity with the analysis. For source-specific analyses, use the stored edge `source` attribute to identify the relationships of interest.

Cite scKNIFE and the biological resources used in the analysis. See the [resource references](../../../docs/graph-resources.md) and [graph selection guide](../../../docs/graphs.md). Resource-specific terms accompany the biological data; the scKNIFE software license applies to the software.
