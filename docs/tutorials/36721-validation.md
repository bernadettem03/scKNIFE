---
orphan: true
---

# Optional: rank sensitivity and H projection

Fit on 29,376 cells and project the other 7,345 onto each learned basis with nonnegative least squares.

## Compare ranks

```python
import scknife

errors = scknife.optional.rank_projection(scknife.tutorial.load_fit(), scknife.tutorial.load_graph(), ranks=(8,16,24))
scknife.optional.plot_projection(errors)
```

```{image} ../assets/tutorial-36721/rank_sensitivity.png
:class: compact-scientific-plot
:align: center
:alt: Held-out reconstruction error at three fitted ranks
```

Relative reconstruction errors are 0.61355, 0.60027 and 0.59354 at ranks 8, 16 and 24, respectively. Because this split overlaps gene selection and tuning, use it as a projection example rather than an independent performance estimate.
