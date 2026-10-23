# Downstream analysis of 10,000 cells

Complete {doc}`tutorials/10000-setup`, then examine {doc}`tutorials/10000-annotation` before exploring biological activities. These analyses reuse the same fitted model and node-by-cell activity scores. Source `cellType` labels provide grouping metadata; scKNIFE does not assign new cell-type labels in this workflow.

1. **Metabolic analysis:** follow {doc}`tutorials/10000-metabolic` for task activity, UMAP patterns, linked reactions and genes, and descriptive On-minus-Pre differences.
2. **Factor interpretation:** follow {doc}`tutorials/10000-factors` to summarize the 16 factors using four consensus clusters and representative biological node groups.
3. **Cell–cell signaling:** follow {doc}`tutorials/10000-signaling` for sender–receiver scores, VEGF examples, permutation FDR, centrality and descriptive condition contrasts.

Run each page's code in order, or use its matching script from {doc}`quickstart`. Each page places code, generated figures and explanations together. For additional questions, {doc}`tutorials/10000-states` illustrates cell-state program summaries. An analysis requiring entities absent from the fitted graph requires an appropriate graph and a new fit.
