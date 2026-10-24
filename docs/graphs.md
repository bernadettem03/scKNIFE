# Graph selection and unification

Choose the graph before tuning and keep that exact graph fixed for the final fit.

## Three ways to choose a graph

1. **Use a supplied unified graph.** Choose Version 1, used throughout the main tutorial, or Version 2, which extends Version 1 with Human-GEM relationships.
2. **Use an individual graph collection.** Choose from the seven preserved collections in `graphs/individual/`: Hallmark, CellChat, GO, markers, cell states, reactions, and metabolic tasks. These are focused collections; they can include relationships from multiple upstream resources, as listed below.
3. **Build a custom graph.** Use the graph-building APIs with source files you obtain for your analysis. Select the biological relationships relevant to your question, use consistent gene identifiers, and record the resulting graph checksum. See {doc}`api/graphs`. A custom graph need not reproduce either supplied unified graph.

### Source terms for supplied and custom graphs

The scKNIFE software license applies to the software. Biological relationships retain their applicable upstream terms whether distributed as individual graphs or as Version 1 or Version 2. Offering a custom-build option does not replace the redistribution requirements for graph files supplied in this repository.

When users obtain upstream inputs themselves and build a graph locally, the input distribution occurs through those upstream sources; users still need to follow their applicable terms. See {doc}`graph-resources` for citations, recorded terms and incomplete source records. Availability of a file here is not a claim that every upstream resource has unrestricted reuse terms.

## Supplied unified graphs

| Graph | Intended use | Biological programs | Nodes | Edges |
|---|---|---:|---:|---:|
| **Version 1** | Version 1: tutorial graph | 7,521 | 22,951 | 216,079 |
| **Version 2** | Version 1 plus Human-GEM | 7,668 | 39,247 | 336,225 |

The biological-program count is a sum of category memberships, not the number of graph nodes. Both versions contain 50 Hallmark sets, 218 metabolic tasks, 787 reaction-associated programs, 95 cell-type marker sets, 83 cell-state programs, 6,036 GO biological processes, and 252 CellChat pathways. Version 2 additionally contains 147 Human-GEM pathways. It does not add KEGG or Cytopus collections.

```python
import networkx as nx
G_v1 = nx.read_graphml("graphs/version_1/unified_graph.graphml")
G_v2 = nx.read_graphml("graphs/version_2/unified_graph.graphml")
```

The exact composition, structural counts, and SHA-256 checksums are recorded in `graphs/manifest.json`. Version 1 is intentionally frozen. Version 2 was constructed by merging the Human-GEM component used in the manuscript analysis into that frozen Version 1 graph; this prevents later changes in upstream marker resources from altering the analysis graph.

## Focused graph collections

`graphs/individual/` provides preserved GraphML collections focused on Hallmark sets, CellChat, GO, markers, cell states, reactions, and metabolic tasks. Each collection retains relationships among its included nodes from the source labels shown below. The collection name denotes its focus; the source inventory describes the full contents.

| Collection | Nodes | Edges | Recorded relationship sources |
|---|---:|---:|---|
| `cell_state` | 498 | 696 | `cellchat`, `cellmarker2.0` |
| `cellchat` | 1,160 | 4,144 | `ARCHS4 TFs Coexp Edges`, `breast-PPI`, `cellchat` |
| `enrichr_go` | 5,824 | 56,870 | `ARCHS4 TFs Coexp Edges`, `G_process`, `Human-GEM`, `breast-PPI`, `cellchat`, `rxn_gene` |
| `hallmark` | 4,434 | 50,474 | `ARCHS4 TFs Coexp Edges`, `Human-GEM`, `breast-PPI`, `cellchat`, `hallmark`, `rxn_gene` |
| `markers` | 1,064 | 3,745 | `ARCHS4 TFs Coexp Edges`, `breast-PPI`, `cellchat`, `marker_graph` |
| `rxn_by_gene` | 1,705 | 5,992 | `ARCHS4 TFs Coexp Edges`, `Human-GEM`, `breast-PPI`, `rxn_gene` |
| `task_by_rxn` | 1,148 | 10,989 | `ARCHS4 TFs Coexp Edges`, `Human-GEM`, `KEGG 2021 Human Pathway Gene Edges`, `breast-PPI`, `metabolic_task`, `rxn_gene` |

Each collection's README records its SHA-256 checksum and edge-source membership counts. Multi-source edges contribute to each applicable source category. These focused collections have their own file identities and source composition; select the desired collection explicitly for each analysis.

```python
import networkx as nx
G = nx.read_graphml("graphs/individual/hallmark/hallmark.graphml")
```

Keep the selected graph fixed during tuning and final fitting, and report its file name and checksum. Expression genes and graph nodes should use matching identifiers. See {doc}`graph-resources` for references to the underlying resources.

## Rebuilding and resource provenance

Install `.[graphs]` to rebuild marker graphs; the included frozen GraphML files do not require decoupler or network downloads. The optional builder uses the decoupler 1.x resource API. Live resource retrieval can change over time and is not equivalent to the frozen benchmarking graph.

For reproducible reuse, record the graph checksum together with the source release/accession, retrieval information, citation, and resource-specific terms. The scKNIFE software and biological resources are described separately in {doc}`citation` and {doc}`graph-resources`.
