# GO-focused graph collection

This preserved GraphML collection provides a focused view of the biological network. Its name identifies the collection focus; the file retains relationships from the source labels listed below.

## File identity

- File: `enrichr_go.graphml`
- Nodes: 5,824
- Edges: 56,870
- SHA-256: `6f141905a5bc98aa21e3129caf0edd93aa040da43ff860d3222430b67ad2017f`

## Relationship sources

Counts below are edge-source memberships. An edge with multiple source labels contributes once to each label, so these counts can sum to more than the number of edges.

| Recorded source label | Edge memberships | Resource reference |
|---|---:|---|
| `ARCHS4 TFs Coexp Edges` | 17,065 | [ARCHS4](https://doi.org/10.1038/s41467-018-03751-6) |
| `G_process` | 9,821 | [Enrichr GO collection](https://doi.org/10.1093/nar/gkw377) |
| `Human-GEM` | 37 | [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) |
| `breast-PPI` | 29,675 | Source identifier retained from the supplied graph |
| `cellchat` | 687 | [CellChat](https://doi.org/10.1038/s41467-021-21246-9) |
| `rxn_gene` | 23 | [Metabolic task/reaction resources](https://github.com/earmingol/scCellFie) |

## Use and citation

Load the complete collection with `networkx.read_graphml()`. Select this graph before tuning and keep its node order and file identity with the analysis. For source-specific analyses, use the stored edge `source` attribute to identify the relationships of interest.

Cite scKNIFE and the biological resources used in the analysis. See the [resource references](../../../docs/graph-resources.md) and [graph selection guide](../../../docs/graphs.md). Resource-specific terms accompany the biological data; the scKNIFE software license applies to the software.
