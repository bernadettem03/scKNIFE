# Annotation review: 10,000 cells

Follow the dedicated {doc}`tutorials/10000-annotation` page for cell annotation and marker review, with its code, plots and interpretation together.

The current {doc}`tutorial-10000` documents primary refined annotations, source-reference comparisons, marker evidence and rare-cell sensitivity on the new 10,000-cell fit. Counts and cluster numbers from the smaller development fit are not carried forward.

The annotation rules use measured expression and cluster membership, not source-reference labels. Original graph-score methods are retained as diagnostics in `scripts/run_tutorial_10000_annotation.py`; their agreement is not an accuracy test. Review `breast_annotation/cluster_decisions.csv`, `rare_refinement/reference_metrics.csv`, `rare_refinement/threshold_stability.csv` and `rare_refinement/stromal_cluster_detection.csv` under the new output directory.
