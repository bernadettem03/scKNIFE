import scknife

fit = scknife.tutorial.load_fit()
graph = scknife.tutorial.load_graph()
activity = scknife.tutorial.load_activities(fit)
scores = scknife.activity.node_scores(activity, graph, "GO_term")
score_name = "GO-node activity"

comparison = dict(
    cell_type_key="cellType", condition_key="timepoint",
    subset={"expansion": "E"}, contrast=("Pre", "On"), min_cells=20,
)
go_de = scknife.activity.compare_activity(
    fit, scores, kind="GO_term", score_name=score_name, **comparison
)

cell_type = "Cancer_cell"  # Other options: "T_cell", "Myeloid_cell"
cutoffs = dict(q_cutoff=0.05, min_d=0.3, min_log2fc=0.25)
passing = scknife.activity.select_results(go_de, cell_type=cell_type, **cutoffs)
columns = ["cell_type", "feature", "log2FC", "cohens_d", "adj_p_value"]
passing[columns].head().style.format({"log2FC": "{:.3f}", "cohens_d": "{:.3f}", "adj_p_value": "{:.2e}"})

scknife.activity.plot_effects(go_de, **cutoffs)

scknife.activity.plot_volcano(go_de, cell_type=cell_type, **cutoffs)
