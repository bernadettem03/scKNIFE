# scKNIFE

Fit a graph-guided model to single-cell expression and explore biological activities. The tutorial uses all 36,721 TNBC cells from six patients in the dataset.

```{image} assets/scknife-workflow-v4.png
:alt: Expression and a biological graph are used to fit cell factors and reconstruct node activities.
:align: center
```

## Getting started

[Install scKNIFE](installation.md), [prepare the inputs and fit the model](tutorials/36721-setup.md), then [cluster the cells](tutorials/36721-clustering.md).

## Analysis tasks

Choose any task using the saved model:

1. [Factor interpretation](tutorials/36721-factors.md) — group factors and inspect their biological programs.
2. [Metabolic tasks](tutorials/36721-metabolic.md) — explore MEAI and compare conditions.
3. [GO analysis](tutorials/36721-go.md) — compare biological-process activities.
4. [Cell–cell signaling](tutorials/36721-signaling.md) — explore sender–receiver scores and ligand–receptor pairs.

```{toctree}
:hidden:
:caption: Getting started
:maxdepth: 1

Installation <installation>
Run scKNIFE <tutorials/36721-setup>
Clustering <tutorials/36721-clustering>
```

```{toctree}
:hidden:
:caption: Analysis tasks
:maxdepth: 1

Factor interpretation <tutorials/36721-factors>
Metabolic tasks <tutorials/36721-metabolic>
GO analysis <tutorials/36721-go>
Cell–cell signaling <tutorials/36721-signaling>
```

## More examples

[Cell-state programs](tutorials/36721-states.md) · [Rank sensitivity and H projection](tutorials/36721-validation.md)

{download}`Download the notebook with results <downloads/scKNIFE_36721_tutorial.ipynb>` · [Graph guide](graphs.md) · [API reference](api/index.md) · [Resources](graph-resources.md) · [Citation](citation.md) · [Troubleshooting](troubleshooting.md)
