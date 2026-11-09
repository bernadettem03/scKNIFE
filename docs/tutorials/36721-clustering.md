# Clustering

Cluster cells using their fitted H values. For the displayed UMAP, we choose 6 principal components, 15 neighbors and Leiden resolution 0.5. In `run_h_graph_clustering_fixed`, “fixed” means one configuration per call. You can change any of these settings; they are separate from the model hyperparameters selected during fitting.

```python
import scknife

fit = scknife.tutorial.load_fit()
fit, clusters = scknife.run_h_graph_clustering_fixed(
    fit, pca_dims=6, n_neighbors=15, resolution=0.5, random_state=42)
clusters.to_csv(scknife.tutorial.output_dir() / "clustering.csv", index=False)
scknife.tutorial.save_fit(fit)
scknife.tutorial.plot_clusters(fit)
```

```{image} ../assets/tutorial-36721/cellular_organization/organization_umap.png
:class: compact-scientific-plot
:align: center
:alt: Each point is one of the 36,721 cells
```

Each point is one of the 36,721 cells. The 19 colors identify the Leiden clusters in the fitted factor representation at resolution 0.5.

## Compare clustering settings

To try several settings, use the grid-search function on a copy of the fit. It returns the results for every combination and selects the one with the highest silhouette score. To use a different combination, enter its values in the single-configuration call above.

```python
# trial, grid, best = scknife.run_h_graph_clustering(
#     fit.copy(), pca_dims=[6], neighbors=[15], resolutions=[0.3, 0.5, 0.8],
#     random_state=42)
# grid[["pca_dims", "n_neighbors", "resolution", "n_clusters", "silhouette"]]
```
