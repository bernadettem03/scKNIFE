---
orphan: true
---

# Graph selection

Choose a graph before tuning and use it for the final fit. Expression genes and graph nodes should use matching gene identifiers.

## Choose a graph

- **Version 1:** the unified graph used in this tutorial.
- **Version 2:** Version 1 extended with Human-GEM relationships.
- **Individual graphs:** seven collections for focused analyses.
- **Custom graphs:** combine relationships relevant to your analysis using the {doc}`api/graphs` tools.

## Unified graphs

| Graph | Contents | Nodes | Edges |
|---|---|---:|---:|
| **Version 1** | Hallmark, metabolic tasks, reactions, cell-type markers, cell states, GO and CellChat | 22,951 | 216,079 |
| **Version 2** | Version 1 plus Human-GEM | 39,247 | 336,225 |

```python
import networkx as nx
G_v1 = nx.read_graphml("graphs/version_1/unified_graph.graphml")
G_v2 = nx.read_graphml("graphs/version_2/unified_graph.graphml")
```

## Individual graphs

Choose a collection from `graphs/individual/`. Each includes the relationships connecting its selected nodes.

| Folder | Focus |
|---|---|
| `cell_state` | Cell states |
| `cellchat` | CellChat signaling |
| `enrichr_go` | GO biological processes |
| `hallmark` | Hallmark gene sets |
| `markers` | Cell-type markers |
| `rxn_by_gene` | Reactions and genes |
| `task_by_rxn` | Metabolic tasks and reactions |

For example, load the Hallmark graph:

```python
import networkx as nx
G = nx.read_graphml("graphs/individual/hallmark/hallmark.graphml")
```

The [graph files and inventories](https://github.com/bernadettem03/scKNIFE/tree/main/graphs) provide details for each collection. See {doc}`graph-resources` for citations and license notices.

## Build a custom graph

Use the {doc}`api/graphs` tools with your own source files. To use the marker-resource builder, install `.[graphs]`; this extra is not needed to load the supplied GraphML files.
