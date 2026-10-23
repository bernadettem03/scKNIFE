# scKNIFE

Graph-regularized non-negative matrix factorization for interpretable single-cell analysis.

Connect expression data with biological knowledge graphs to learn cell factors, reconstruct node activity, and explore cellular programs.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Get started
:link: quickstart
:link-type: doc
Install scKNIFE, prepare the source files and run the 10,000-cell workflow.
:::
:::{grid-item-card} Run the tutorial
:link: tutorial-10000
:link-type: doc
Work through 10,000 cells: fitting, retention of cell-type organization, metabolism and signaling, with source-reference comparisons.
:::
:::{grid-item-card} Choose a graph
:link: graphs
:link-type: doc
Explore the two unified graphs and seven focused collections, with source inventories.
:::
:::{grid-item-card} Explore the API
:link: api/index
:link-type: doc
Find functions for fitting, clustering, annotation, projection, and signaling.
:::
::::

## From expression to biological activity

```{image} assets/scknife-workflow-v3.png
:alt: scKNIFE analysis workflow connecting expression, knowledge graphs, factorization and downstream analysis
:align: center
```

1. **Install and prepare inputs.** Obtain expression, metadata and the appropriate graph.
2. **Run full scKNIFE.** Sample 10,000 cells, fit W and H, and reconstruct node activities.
3. **Check cell-type retention.** Compare learned clusters with known cell types to assess cell-type retention.
4. **Explore downstream tasks.** Analyze metabolism, interpret factors, then investigate cell–cell signaling.

## A reproducible starting point

The {doc}`tutorial-10000` is the main end-to-end walkthrough. It uses one stratified sample with authentic donor, timepoint and cell-type metadata, and explains the executable code and figures at every stage. It requires the external cohort inputs described in the guide.

Follow {doc}`tutorial-inputs` to obtain and prepare the source data, then use {doc}`quickstart` to generate every result. The quick start and detailed tutorial run the same 10,000-cell analysis.

See {doc}`installation` for the Python environment and {doc}`graph-resources` for biological-resource citations.

```{toctree}
:hidden:
:caption: 10,000-cell tutorial
:maxdepth: 1

Overview <tutorial-10000>
1. Installation and inputs <installation>
2. Run full scKNIFE <tutorials/10000-setup>
3. Cell-type retention <tutorials/10000-annotation>
4. Metabolic analysis <tutorials/10000-metabolic>
5. Factor interpretation <tutorials/10000-factors>
6. Cell–cell signaling <tutorials/10000-signaling>
7. Robustness and exports <tutorials/10000-validation>
```

```{toctree}
:hidden:
:caption: Optional analyses and diagnostics
:maxdepth: 1

Other activity-based analyses <tutorials/10000-states>
Marker-based reassignment <tutorials/10000-diagnostic>
```

```{toctree}
:hidden:
:caption: Reference and resources
:maxdepth: 1

tutorial-inputs
quickstart
concepts
graphs
api/index
graph-resources
citation
release-checklist
troubleshooting
```

```{toctree}
:hidden:
:caption: Additional reference
:maxdepth: 1

Tutorial guide <tutorials/index>
Annotation diagnostics <annotation-audit>
tuning
pipeline
projection
downstream
signaling
```
