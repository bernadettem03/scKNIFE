# Quick start

Complete {doc}`installation` and {doc}`tutorial-inputs`. Activate the environment and set `SCKNIFE_COUNTS`, `SCKNIFE_METADATA` and `SCKNIFE_INTERACTION_DB`.

```bash
export SCKNIFE_TUTORIAL_OUTPUT="$PWD/tutorial/outputs_10000"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export NUMBA_NUM_THREADS=2
python scripts/run_tutorial_10000.py
python scripts/run_tutorial_10000_organization.py
python scripts/run_tutorial_10000_metabolic.py
python scripts/run_tutorial_10000_go.py
python scripts/run_tutorial_10000_factors.py
python scripts/run_tutorial_10000_signaling.py
```

The first command samples 10,000 cells, tunes parameters on a representative 2,000-cell subset, refits all 10,000 cells with the selected settings and reconstructs activities. The following commands assess cell-type retention, test metabolic-task and GO activity, interpret factors at k=4, and analyze signaling. Each command matches the code on its corresponding tutorial page.

The original `cellType` metadata provides grouping and comparison labels. Cell-state activity analysis is an optional extension; it is not required by these commands.

Use a new output directory for a fresh run. Signaling performs 1,000 permutations and resumes only compatible checkpoints. See {doc}`tutorial-10000` for the figures and explanations.

## Optional: rank sensitivity and projection

After the main analyses, run this example to project held-out cells onto a fixed learned basis:

```bash
python scripts/run_tutorial_10000_validation.py
```

See {doc}`tutorials/10000-validation` for the code and plot.
