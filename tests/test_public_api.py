def test_public_api_imports():
    from scknife import (run_full_pipeline, run_h_graph_clustering,
        run_h_graph_clustering_fixed,
        run_celltype_assignment, run_cellstate_assignment, run_activity_lr)
    assert all(callable(x) for x in [run_full_pipeline, run_h_graph_clustering,
        run_h_graph_clustering_fixed, run_celltype_assignment,
        run_cellstate_assignment, run_activity_lr])


def test_assignment_helpers_on_synthetic_data():
    import pandas as pd
    from scknife.celltype import compute_cluster_marker_means, assign_cell_types
    from scknife.cellstate import build_state_mapping_table

    activity = pd.DataFrame(
        [[4.0, 3.0, 0.1, 0.2], [0.2, 0.1, 4.0, 3.0]],
        index=["T cell", "B cell"], columns=["c1", "c2", "c3", "c4"]
    )
    labels = pd.Series([0, 0, 1, 1], index=activity.columns)
    means = compute_cluster_marker_means(activity, activity.index.tolist(), labels)
    assigned, scores = assign_cell_types(means)
    assert assigned.to_dict() == {0: "T cell", 1: "B cell"}
    assert (scores > 0).all()

    table = build_state_mapping_table({
        0: {"cell_state": "active", "score": 2.0, "n_cells": 2},
        1: {"cell_state": "resting", "score": 1.0, "n_cells": 2},
    })
    assert table["Assigned_CellState"].tolist() == ["active", "resting"]
