import scknife

fit = scknife.tutorial.load_fit()
activity_path = scknife.tutorial.output_dir() / "activity_scores.parquet"
signals = scknife.communication.score_signaling(fit, activity_path=activity_path, groupby="cellType")

scknife.communication.plot_matrix(signals)

scknife.communication.plot_strength(signals)

scknife.communication.plot_contributions(signals)

fdr = scknife.communication.permutation_fdr(signals, n_permutations=1000)
scknife.communication.plot_fdr(fdr)

scknife.communication.plot_network(signals)

central = scknife.communication.centrality(signals)
scknife.communication.plot_centrality(central)

pathway = "VEGF"
# Other pathways: "CXCL", "TGFb", "WNT"
matrix, pairs = scknife.communication.pathway_scores(signals, pathway)
scknife.communication.plot_pathway(matrix, pathway)

scknife.communication.plot_pathway_pairs(pairs, pathway)
