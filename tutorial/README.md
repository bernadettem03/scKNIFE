# Complete 10,000-cell tutorial

Open [the complete notebook](scKNIFE_Bassez_Cohort_B_complete_tutorial.ipynb) or follow [the step-by-step documentation](../docs/tutorial-10000.md). Both use the same representative 10,000-cell sample and generate their figures from the fitted model.

## Prepare and run

Install Python 3.9.21 and `pip install -e '.[tutorial]'` from the repository root. Follow [input preparation](../docs/tutorial-inputs.md) for source counts, metadata and the directional interaction table. Set `SCKNIFE_COUNTS`, `SCKNIFE_METADATA` and `SCKNIFE_INTERACTION_DB`, then run:

```bash
python scripts/run_tutorial.py
```

The executed notebook and results are saved under `tutorial/outputs_10000/`, or the directory selected by `SCKNIFE_TUTORIAL_OUTPUT`. Alternatively, use the six commands in [Quick start](../docs/quickstart.md).

## Workflow

Fit scKNIFE and reconstruct activities → examine cell-type retention → metabolic analysis → factor interpretation → cell–cell signaling → robustness and exports. The primary fit uses 16 factors; interpretation groups factors into four consensus clusters. Source cell-type metadata remain attached to the cells. VEGF and On minus Pre are worked signaling examples; E/NE denotes clonotype expansion/non-expansion.

## Separate regression example

[The 1,000-cell regression notebook](scKNIFE_regression_1000.ipynb) uses the bundled synthetic test fixture and includes synthetic signaling checks. It is separate from the main biological tutorial:

```bash
python scripts/run_tutorial.py --regression
```

Regression outputs default to `tutorial/outputs/`. Use a separate `SCKNIFE_TUTORIAL_OUTPUT` if that variable is already set. See [fixture provenance](data/README.md) for the bundled data records.
