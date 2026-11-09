# Cell–cell signaling

Explore sender–receiver scores, ligand–receptor contributions and a pathway of interest. The first run prepares the interaction table using R and internet access.

```python
import scknife

fit = scknife.tutorial.load_fit()
activity_path = scknife.tutorial.output_dir() / "activity_scores.parquet"
signals = scknife.communication.score_signaling(fit, activity_path=activity_path, groupby="cellType")
```

## Sender–receiver activity

```python
scknife.communication.plot_matrix(signals)
```

```{image} ../assets/tutorial-36721/signaling/signaling_heatmap.png
:class: compact-scientific-plot signaling-single
:align: center
:alt: Rows send and columns receive
```

Rows send and columns receive. Color shows the pooled activity-weighted interaction score, including self-links.

## Outgoing and incoming strength

```python
scknife.communication.plot_strength(signals)
```

```{image} ../assets/tutorial-36721/signaling/general/01_sender_receiver.png
:class: compact-scientific-plot signaling-wide
:align: center
:alt: Outgoing and incoming scores summarize a cell type’s activity as a sender and receiver
```

Outgoing and incoming scores summarize a cell type’s activity as a sender and receiver.

## Ligand–receptor contributions

```python
scknife.communication.plot_contributions(signals)
```

```{image} ../assets/tutorial-36721/signaling/signaling_lr_bubbles.png
:class: compact-scientific-plot signaling-wide
:align: center
:alt: Bubble size shows each component’s weighted contribution within the selected sender–receiver pairs
```

Bubble size shows each component’s weighted contribution within the selected sender–receiver pairs.

## Permutation support

```python
fdr = scknife.communication.permutation_fdr(signals, n_permutations=1000)
scknife.communication.plot_fdr(fdr)
```

```{image} ../assets/tutorial-36721/signaling/inference/permutation_fdr.png
:class: compact-scientific-plot signaling-single
:align: center
:alt: Color shows −log₁₀ of the BH-adjusted p-value across sender–receiver totals
```

Color shows −log₁₀ of the BH-adjusted p-value across sender–receiver totals. In this cohort, 25 of the 64 pairs have adjusted p < 0.05. Shuffling cell labels retains group sizes and uses the observed specificity.

## Signaling network

```python
scknife.communication.plot_network(signals)
```

```{image} ../assets/tutorial-36721/signaling/general/06_hierarchical_network.png
:class: compact-scientific-plot signaling-single
:align: center
:alt: Arrows show signaling direction for up to 40 positive links
```

Arrows show signaling direction for up to 40 positive links. Wider lines indicate stronger scores.

## Network centrality

```python
central = scknife.communication.centrality(signals)
scknife.communication.plot_centrality(central)
```

```{image} ../assets/tutorial-36721/signaling/inference/network_centrality.png
:class: compact-scientific-plot signaling-wide
:align: center
:alt: Strength, betweenness and PageRank summarize the network after self-links are removed
```

Strength, betweenness and PageRank summarize the network after self-links are removed.

## Pathway example: VEGF

```python
pathway = "VEGF"
# Other pathways: "CXCL", "TGFb", "WNT"
matrix, pairs = scknife.communication.pathway_scores(signals, pathway)
scknife.communication.plot_pathway(matrix, pathway)
```

```{image} ../assets/tutorial-36721/signaling/pathways/VEGF/03_communication_heatmap.png
:class: compact-scientific-plot signaling-single
:align: center
:alt: Compare sender–receiver scores within VEGF
```

Compare sender–receiver scores within VEGF. The permutation test above applies to global totals, not to individual pathways.

## Pathway contributors

```python
scknife.communication.plot_pathway_pairs(pairs, pathway)
```

```{image} ../assets/tutorial-36721/signaling/pathways/VEGF/05_top_lr_pairs.png
:class: compact-scientific-plot signaling-single
:align: center
:alt: Bars rank ligand–receptor pairs by total weighted score in the selected pathway
```

Bars rank ligand–receptor pairs by total weighted score in the selected pathway.

{download}`Download this example <../downloads/signaling_walkthrough.py>`.
