# Hyperparameter tuning before fitting

The [main setup page](tutorials/10000-setup.md) tunes model rank, graph penalty and W/H penalties before fitting all 10,000 cells. It contains the executable code, subset-composition plot, parameter-selection plot and measured results.

The representative tuning subset contains 2,000 cells: 1,500 fitting cells and 500 validation cells, sampled across patient, timepoint and source cell type. Genes are selected on fitting cells only. Selection uses reconstruction of validation expression with a fixed learned basis, averaged across two initialization seeds. Source labels are not optimization targets.

The search compares ranks 8, 16 and 24 jointly with graph penalties 0.1, 1 and 10, then varies W and H penalties sequentially. The smallest rank within 1% of the best initial-grid error is retained. The selected parameters and gene list are saved before the full refit. These choices are specific to this data, graph, preprocessing and search budget.

The later {doc}`tutorials/10000-validation` page remains an optional projection demonstration. It reuses the selected penalties and is not an independent test of the tuning choices.
