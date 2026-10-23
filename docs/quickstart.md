# Quick start

Complete {doc}`installation` and {doc}`tutorial-inputs`. Activate the environment and set `SCKNIFE_COUNTS`, `SCKNIFE_METADATA` and `SCKNIFE_INTERACTION_DB`. Run from the repository root.

```bash
export SCKNIFE_TUTORIAL_OUTPUT="$PWD/tutorial/outputs_10000"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export NUMBA_NUM_THREADS=2
python scripts/run_tutorial_10000.py
python scripts/run_tutorial_10000_organization.py
python scripts/run_tutorial_10000_metabolic.py
python scripts/run_tutorial_10000_factors.py
python scripts/run_tutorial_10000_signaling.py
python scripts/run_tutorial_10000_validation.py
```

The first command samples 10,000 cells, fits the model and reconstructs activities. The following commands assess cell-type retention, explore metabolic activity, interpret factors at k=4, analyze signaling and run the rank/projection diagnostic. Each command matches the code on its corresponding tutorial page.

The original `cellType` metadata provides grouping and comparison labels. Marker-based reassignment and cell-state analysis are optional extensions; neither is required by these commands.

Use a new output directory for a fresh run. Signaling performs 1,000 permutations and resumes only compatible checkpoints. See {doc}`tutorial-10000` for the figures and explanations.
