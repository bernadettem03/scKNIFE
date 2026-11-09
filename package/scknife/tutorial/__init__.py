"""Convenience functions for the scKNIFE tutorial."""

from .core import (
    output_dir,
    load_graph,
    build_graph,
    marker_graph,
    optimize_hyperparameters,
    fit_and_score,
    load_fit,
    load_activities,
    save_fit,
    figure,
    sample_bassez,
    plot_sample,
    prepare_tuning,
    select_parameters,
    plot_convergence,
    plot_clusters,
)

from .activity import plot_task_umap

__all__ = [
    "output_dir",
    "load_graph",
    "build_graph",
    "marker_graph",
    "optimize_hyperparameters",
    "fit_and_score",
    "load_fit",
    "load_activities",
    "save_fit",
    "figure",
    "sample_bassez",
    "plot_sample",
    "prepare_tuning",
    "select_parameters",
    "plot_convergence",
    "plot_clusters",
    "plot_task_umap",
]
