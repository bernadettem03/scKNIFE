# Held-out projection within 10,000 cells

Follow the dedicated {doc}`tutorials/10000-validation` page for held-out projection, with its code, plots and interpretation together.

The executed diagnostic in {doc}`tutorial-10000` partitions the selected 10,000 cells into 8,000 training cells and 2,000 disjoint held-out cells using seed 42. Each rank is fitted on the training cells; nonnegative least squares estimates H for held-out cells while keeping the learned expression basis `Rg @ W` fixed.

```bash
python scripts/run_tutorial_10000_validation.py
```

The script reloads `tutorial/outputs_10000/fitted_subset.h5ad`, fits ranks 8, 16 and 24, and saves `rank_sensitivity.csv`, `heldout_projected_factors.csv`, `rank_sensitivity.png` and `validation_report.json`. The displayed main UMAP and biological summaries come from the separate final model fitted on all 10,000 cells.

Gene selection and parameter tuning used a subset that overlaps this exploratory split, so the diagnostic is exploratory rather than an unbiased generalization benchmark. Keep graph-node order, gene order and preprocessing consistent. No additional cohort cells are projected here. The package's general `h_projection` API is described in {doc}`api/index`.
