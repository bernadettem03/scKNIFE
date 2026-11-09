import scknife

sample = scknife.tutorial.sample_bassez(n_cells=None, seed=42)

# graph = scknife.tutorial.build_graph(["hallmark", "enrichr_go", "rxn_by_gene", "task_by_rxn"])
graph = scknife.tutorial.load_graph("graphs/version_1/unified_graph.graphml")
# graph = scknife.tutorial.load_graph("graphs/version_2/unified_graph.graphml")

adata, params = scknife.tutorial.optimize_hyperparameters(
    sample, graph, ranks=(8, 16, 24), subset=2000, validation=500, seed=42
)
params

fit = scknife.tutorial.fit_and_score(adata, graph, params)
