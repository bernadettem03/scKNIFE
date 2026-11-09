import scknife

fit = scknife.tutorial.load_fit()
graph = scknife.tutorial.load_graph()
activity = scknife.tutorial.load_activities(fit)
scores = scknife.activity.calculate_meai(activity, graph)
score_name = "MEAI"
# scores, score_name = scknife.activity.node_scores(activity, graph, "task"), "Task-node activity"

comparison = dict(
    cell_type_key="cellType", condition_key="timepoint",
    subset={"expansion": "E"}, contrast=("Pre", "On"), min_cells=20,
)
task_de = scknife.activity.compare_activity(
    fit, scores, kind="task", score_name=score_name, **comparison
)

cell_type = "Cancer_cell"  # Other options: "T_cell", "Myeloid_cell"
cutoffs = dict(q_cutoff=0.05, min_d=0.3, min_log2fc=0.25)
passing = scknife.activity.select_results(task_de, cell_type=cell_type, **cutoffs)
columns = ["cell_type", "feature", "log2FC", "cohens_d", "adj_p_value"]
passing[columns].head().style.format({"log2FC": "{:.3f}", "cohens_d": "{:.3f}", "adj_p_value": "{:.2e}"})

scknife.activity.plot_effects(task_de, **cutoffs)

scknife.activity.plot_volcano(task_de, cell_type=cell_type, **cutoffs)

tasks = ["Glucose to lactate conversion", "GMP salvage from guanine"]
# For one panel: tasks = "Mevalonate synthesis"
scknife.tutorial.plot_task_umap(fit, scores, task=tasks, score_name=score_name)

scknife.activity.plot_task_violin(fit, scores, task="GMP salvage from guanine",
                     cell_type=cell_type, score_name=score_name, **comparison)
