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

| Collection | Relationships it contributes |
|---|---|
| `hallmark` | Genes linked to MSigDB Hallmark programs for factor interpretation. |
| `enrichr_go` | Genes linked to functional process sets. Check the set names when choosing a process to study. |
| `rxn_by_gene` | Gene–reaction associations for interpreting metabolic activity. |
| `task_by_rxn` | Metabolic task-support relationships. Pair with `rxn_by_gene` for reaction context. |
| `cellchat` | Ligand and receptor components linked through signaling interactions. |
| `markers` | Marker genes linked to cell types, including breast epithelial and tumor-microenvironment populations. |
| `cell_state` | Genes linked to state or population programs from the breast-tissue CellMarker panel. |

Each GraphML is a focused collection: it also retains other relationships between the selected nodes. For example, merging a different marker panel does not replace breast-tissue protein interactions present in other collections. For GO differential analysis, the fitted graph must contain nodes with `kind="GO_term"`; the individual process collection alone is not a substitute for all GO terms in Version 1.

For example, load the Hallmark graph:

```python
import networkx as nx
G = nx.read_graphml("graphs/individual/hallmark/hallmark.graphml")
```

The [graph files and inventories](https://github.com/bernadettem03/scKNIFE/tree/main/graphs) provide details for each collection. See {doc}`graph-resources` for citations and license notices.

## Build a custom graph

Combine the supplied collections you need:

```python
import scknife

graph = scknife.tutorial.build_graph(["hallmark", "enrichr_go", "rxn_by_gene", "task_by_rxn"])
```

For your own prepared GraphML collections, pass their graphs to the same merger:

```python
import scknife

graph = scknife.construct_unified_graph({"pathways": scknife.tutorial.load_graph("my_pathways.graphml"),
                                 "signaling": scknife.tutorial.load_graph("my_signaling.graphml")})
```

Match gene identifiers to the counts and retain the node kinds needed by your analysis. See {doc}`api/graphs` for individual builders. To use the marker-resource builder, install `.[graphs]`; this extra is not needed to load the supplied GraphML files.

(use-markers-for-another-cancer)=
## Marker set

For this dataset, `markers` includes breast epithelial populations and broader immune and stromal populations. The separate `cell_state` collection includes a breast-tissue panel from CellMarker. Select the two inputs independently when constructing a graph for another tissue.

1. In [CellMarker 2.0](https://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/index.html), select the species, tissue and cancer type that match your cohort. For example, use **human → lung → the relevant lung cancer subtype** where records are available. Keep the cell populations and supporting studies relevant to your samples. A curated panel from a tumor-specific atlas is another option. [CellMarker 2.0 publication](https://doi.org/10.1093/nar/gkac947)
2. Save the filtered pairs in a CSV with columns `gene` and `cell_type`: one gene symbol and one cell-type label per row. Use gene identifiers matching the expression matrix. A marker can appear on more than one row when it supports several cell types.
3. Add that panel to your selected collections **before tuning**. Start from individual collections and leave out `markers` to replace the breast marker input. Leave out `cell_state` too unless its programs are appropriate for your cohort.

```python
import scknife

base = scknife.tutorial.build_graph(["hallmark", "rxn_by_gene", "task_by_rxn"])
panel = scknife.tutorial.marker_graph("lung_markers.csv", source="CellMarker 2.0")
graph = scknife.construct_unified_graph({"markers": panel, "metabolism_and_programs": base})
```

Pass `graph` to hyperparameter optimization and the full fit. This changes the biological relationships available to the model; it does not assign cell types to the cells.
