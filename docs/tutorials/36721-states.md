---
orphan: true
---

# Optional: cell-state programs

## Activity across clusters

```python
import scknife

fit = scknife.tutorial.load_fit()
means = scknife.optional.state_programs(fit, scknife.tutorial.load_activities(fit), scknife.tutorial.load_graph(), n=8)
scknife.optional.plot_states(means)
```

```{image} ../assets/tutorial-36721/state_activity_by_cluster.png
:class: compact-scientific-plot
:align: center
:alt: Rows are the eight state programs with highest pooled mean activity
```

Rows are the eight state programs with highest pooled mean activity. Color shows their mean activity in each cluster.
