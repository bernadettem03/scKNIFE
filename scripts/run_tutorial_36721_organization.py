import scknife

fit = scknife.tutorial.load_fit()
fit, clusters = scknife.run_h_graph_clustering_fixed(
    fit, pca_dims=6, n_neighbors=15, resolution=0.5, random_state=42)
clusters.to_csv(scknife.tutorial.output_dir() / "clustering.csv", index=False)
scknife.tutorial.save_fit(fit)
scknife.tutorial.plot_clusters(fit)

# trial, grid, best = scknife.run_h_graph_clustering(
#     fit.copy(), pca_dims=[6], neighbors=[15], resolutions=[0.3, 0.5, 0.8],
#     random_state=42)
# grid[["pca_dims", "n_neighbors", "resolution", "n_clusters", "silhouette"]]
