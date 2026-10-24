# The 10,000-cell scKNIFE tutorial

Fit scKNIFE on 10,000 Bassez cells, explore clustering, then choose an analysis for your question.

## Get started

1. [Install scKNIFE](installation.md) and [prepare the input files](tutorial-inputs.md).
2. [Run scKNIFE](tutorials/10000-setup.md): sample cells, tune parameters, fit the model and save activities.
3. [Explore clustering and cell-type retention](tutorials/10000-annotation.md).

## Choose an analysis

Each task uses the saved model and activities. Choose any task; the analysis pages do not depend on one another.

| Task | What you can do |
| --- | --- |
| [Metabolic tasks](tutorials/10000-metabolic.md) | Calculate MEAI, plot a task and compare conditions. |
| [GO analysis](tutorials/10000-go.md) | Compare GO biological-process activities between conditions. |
| [Factor interpretation](tutorials/10000-factors.md) | Group factors and explore their biological programs. |
| [Cell–cell signaling](tutorials/10000-signaling.md) | Explore sender–receiver scores, ligand–receptor pairs and pathways. |

## More examples

- [Rank sensitivity and H projection](tutorials/10000-validation.md)
- [Cell-state activity](tutorials/10000-states.md)

Other analyses can use activities for entities in your graph. To study additional entities, supply a suitable graph and refit.

## Notebook and scripts

Download the {download}`notebook <downloads/scKNIFE_10000_tutorial.ipynb>` for all examples, or use the {doc}`quickstart` to run individual tasks.
