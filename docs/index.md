# scKNIFE

Graph-regularized non-negative matrix factorization for interpretable single-cell analysis.

Connect expression data with biological knowledge graphs to learn cell factors, reconstruct node activity, and explore cellular programs.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Getting started
:link: quickstart
:link-type: doc
Install scKNIFE, prepare the source files and run the 10,000-cell workflow.
:::
:::{grid-item-card} Run the tutorial
:link: tutorial-10000
:link-type: doc
Fit and cluster 10,000 cells, then choose a biological analysis.
:::
:::{grid-item-card} Choose a graph
:link: graphs
:link-type: doc
Choose Version 1, Version 2 or an individual graph.
:::
:::{grid-item-card} Explore the API
:link: api/index
:link-type: doc
Find functions for fitting, clustering, activity analysis, projection, and signaling.
:::
::::

## From expression to biological activity

```{image} assets/scknife-workflow-v4.png
:alt: scKNIFE analysis workflow connecting expression, knowledge graphs, factorization and downstream analysis
:align: center
```

Start with {doc}`installation`, {doc}`tutorial-inputs` and {doc}`tutorials/10000-setup`. Explore clustering and cell-type retention, then choose metabolic tasks, GO analysis, factor interpretation or signaling from the {doc}`tutorial-10000`.

```{toctree}
:hidden:
:caption: Getting started
:maxdepth: 1

Tutorial overview <tutorial-10000>
Installation <installation>
Input files <tutorial-inputs>
Run scKNIFE <tutorials/10000-setup>
Clustering and cell-type retention <tutorials/10000-annotation>
```

```{toctree}
:hidden:
:caption: Analysis tasks
:maxdepth: 1

Metabolic tasks <tutorials/10000-metabolic>
GO analysis <tutorials/10000-go>
Factor interpretation <tutorials/10000-factors>
Cell–cell signaling <tutorials/10000-signaling>
```

[Quick start](quickstart.md) · [Concepts](concepts.md) · [API reference](api/index.md) · [Resource references](graph-resources.md) · [Citation](citation.md) · [Troubleshooting](troubleshooting.md)
