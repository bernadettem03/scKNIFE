# Synthetic regression fixture

`synthetic_regression_1000.h5ad` contains 1,000 simulated cells and 1,000 genes, with invented counts and four synthetic groups. No patient expression, cell identifiers or metadata were used to generate these data. Gene names are selected from the unchanged Version 1 graph.

Regenerate with `python scripts/generate_synthetic_fixture.py` (seed 42). Counts follow Poisson sampling from random group profiles and cell library multipliers. The CSV contains only synthetic group labels. This is a software test fixture, not a biological benchmark; graph-derived annotations and activities must not be interpreted as biological findings.

The synthetic data are supplied under the repository MIT license. Graph resource terms are recorded separately.

The main [10,000-cell tutorial](../../docs/tutorial-10000.md) continues to use Bassez cohort-2 data obtained from the [original source](https://lambrechtslab.sites.vib.be/en/single-cell). Follow [input preparation](../../docs/tutorial-inputs.md); the real counts and metadata are not bundled here.
