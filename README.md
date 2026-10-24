# scKNIFE

**Connect single-cell gene expression with biological knowledge to reveal interpretable cellular programs.**

[Documentation](https://scknife-docs.readthedocs.io/en/latest/) · [Tutorial notebook](tutorial/scKNIFE_Bassez_Cohort_B_complete_tutorial.ipynb) · [Installation](https://scknife-docs.readthedocs.io/en/latest/installation.html) · [Graph guide](https://scknife-docs.readthedocs.io/en/latest/graphs.html)

[![Tests](https://github.com/bernadettem03/scKNIFE/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/bernadettem03/scKNIFE/actions/workflows/tests.yml)
[![Documentation](https://readthedocs.org/projects/scknife-docs/badge/?version=latest)](https://scknife-docs.readthedocs.io/en/latest/)
[![Software license: MIT](https://img.shields.io/badge/software_license-MIT-blue.svg)](LICENSE)

scKNIFE uses graph-regularized non-negative matrix factorization to learn cell factors and biological-node activities.

<p align="center">
  <img src="docs/assets/scknife-workflow-v4.png" alt="scKNIFE workflow: gene expression and a biological knowledge graph lead to cell factors, node loadings and downstream biological activity analyses" width="850">
</p>

## Getting started

Use Python 3.9.21 and install the tutorial dependencies:

```bash
git clone https://github.com/bernadettem03/scKNIFE.git
cd scKNIFE
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[tutorial]'
```

Prepare the [input files](https://scknife-docs.readthedocs.io/en/latest/tutorial-inputs.html) and set their paths using the [quick start](https://scknife-docs.readthedocs.io/en/latest/quickstart.html). Open the [10,000-cell notebook](tutorial/scKNIFE_Bassez_Cohort_B_complete_tutorial.ipynb) in Jupyter, or run all examples with:

```bash
python scripts/run_tutorial.py
```

The tutorial covers parameter tuning, fitting and [clustering and cell-type retention](https://scknife-docs.readthedocs.io/en/latest/tutorials/10000-annotation.html). Existing cell-type labels are used for interpretation; scKNIFE does not assign new labels. Results are saved in your configured output directory.

## Analysis tasks

Choose any task using the saved fit and activities:

| Task | What you can do |
|---|---|
| [Metabolic tasks](https://scknife-docs.readthedocs.io/en/latest/tutorials/10000-metabolic.html) | Calculate MEAI, plot a task and compare conditions. |
| [GO analysis](https://scknife-docs.readthedocs.io/en/latest/tutorials/10000-go.html) | Compare GO biological-process activities between conditions. |
| [Factor interpretation](https://scknife-docs.readthedocs.io/en/latest/tutorials/10000-factors.html) | Group factors and explore their biological programs. |
| [Cell–cell signaling](https://scknife-docs.readthedocs.io/en/latest/tutorials/10000-signaling.html) | Explore sender–receiver scores, ligand–receptor pairs and pathways. |

See the [tutorial overview](https://scknife-docs.readthedocs.io/en/latest/tutorial-10000.html) for optional examples, including H projection, and the [API reference](https://scknife-docs.readthedocs.io/en/latest/api/index.html) for package functions.

## Choose a graph

| Option | Use |
|---|---|
| **Ready-made Version 1 or Version 2** | Start with a unified biological graph. The tutorial uses Version 1; Version 2 adds Human-GEM. |
| **Seven individual collections** | Select Hallmark, CellChat, GO, markers, cell states, reactions or metabolic-task relationships. |
| **Custom graph** | Build a graph with the supplied APIs and upstream inputs appropriate to your question. |

Explore the [graph guide](https://scknife-docs.readthedocs.io/en/latest/graphs.html) and [supplied graph files](graphs/). Choose the graph before fitting. Analyses that require additional entities need a suitable graph and a new fit; downstream analyses can reuse the resulting activity scores.

## Citation and support

When using scKNIFE, cite the repository and the version or commit used, together with the biological resources included in your analysis. The manuscript citation will be added when available. See [citation guidance](https://scknife-docs.readthedocs.io/en/latest/citation.html) and [resource references](https://scknife-docs.readthedocs.io/en/latest/graph-resources.html).

For questions or reproducible problems, [open a GitHub issue](https://github.com/bernadettem03/scKNIFE/issues). Include your scKNIFE version, environment and a minimal example without private data.

## License

scKNIFE software and documentation use the [MIT License](LICENSE). Supplied biological resources retain their source-specific terms; see [resource references](https://scknife-docs.readthedocs.io/en/latest/graph-resources.html). These terms apply to both individual and unified graphs, and to upstream inputs used for custom graph construction.
