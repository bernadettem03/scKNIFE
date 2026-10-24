---
orphan: true
---

# Quick start

Complete {doc}`installation` and {doc}`tutorial-inputs`. Activate the environment and set `SCKNIFE_COUNTS`, `SCKNIFE_METADATA` and `SCKNIFE_INTERACTION_DB`.

## Run scKNIFE and explore clustering

```bash
export SCKNIFE_TUTORIAL_OUTPUT="$PWD/tutorial/outputs_10000"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export NUMBA_NUM_THREADS=2
python scripts/run_tutorial_10000.py
python scripts/run_tutorial_10000_organization.py
```

The first script samples 10,000 cells, tunes on 2,000 cells, fits the full sample and saves activities. The second displays clusters alongside the existing cell-type labels. Use a new output directory for a fresh fit.

## Analysis tasks

Run any of these tasks using the saved outputs above.

### Metabolic tasks

```bash
python scripts/run_tutorial_10000_metabolic.py
```

See {doc}`tutorials/10000-metabolic` for MEAI, a task violin plot and differential activity.

### GO analysis

```bash
python scripts/run_tutorial_10000_go.py
```

See {doc}`tutorials/10000-go` for GO differential activity.

### Factor interpretation

```bash
python scripts/run_tutorial_10000_factors.py
```

See {doc}`tutorials/10000-factors` for factor groups and biological programs.

### Cell–cell signaling

```bash
python scripts/run_tutorial_10000_signaling.py
```

See {doc}`tutorials/10000-signaling` for communication scores and pathway plots. The 1,000-permutation calculation resumes compatible checkpoints.

## Optional: rank sensitivity and projection

```bash
python scripts/run_tutorial_10000_validation.py
```

See {doc}`tutorials/10000-validation` for projecting held-out cells onto a learned basis.
