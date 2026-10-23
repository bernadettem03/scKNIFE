# Metabolic-task graph collection

This preserved GraphML collection provides a focused view of the biological network. Its name identifies the collection focus; the file retains relationships from the source labels listed below.

## File identity

- File: `task_by_rxn.graphml`
- Nodes: 1,148
- Edges: 10,989
- SHA-256: `6bbcfa958cebd588477f6b5105f47164b4068a92497561b67fc3c1f456ce4e61`

## Relationship sources

Counts below are edge-source memberships. An edge with multiple source labels contributes once to each label, so these counts can sum to more than the number of edges.

| Recorded source label | Edge memberships | Resource reference |
|---|---:|---|
| `ARCHS4 TFs Coexp Edges` | 1,134 | [ARCHS4](https://doi.org/10.1038/s41467-018-03751-6) |
| `Human-GEM` | 6 | [Human-GEM](https://github.com/SysBioChalmers/Human-GEM) |
| `KEGG 2021 Human Pathway Gene Edges` | 17 | [KEGG](https://www.kegg.jp/) |
| `breast-PPI` | 2,678 | Source identifier retained from the supplied graph |
| `metabolic_task` | 7,293 | [scCellFie task resources](https://github.com/earmingol/scCellFie) |
| `rxn_gene` | 21 | [Metabolic task/reaction resources](https://github.com/earmingol/scCellFie) |

## Use and citation

Load the complete collection with `networkx.read_graphml()`. Select this graph before tuning and keep its node order and file identity with the analysis. For source-specific analyses, use the stored edge `source` attribute to identify the relationships of interest.

Cite scKNIFE and the biological resources used in the analysis. See the [resource references](../../../docs/graph-resources.md) and [graph selection guide](../../../docs/graphs.md). Resource-specific terms accompany the biological data; the scKNIFE software license applies to the software.
