# The 10,000-cell scKNIFE tutorial

**Expression + knowledge graph → tune parameters → run scKNIFE → check cell-type retention → explore biological activities.**

Run the model on one representative sample of 10,000 cells, then use its learned factors and reconstructed node activities for downstream questions. scKNIFE does not assign cell-type labels in this workflow. Comparing the learned representation with known cell types assesses whether it retains recognizable biological organization.

## Follow the workflow

1. [Installation and inputs](installation.md): install the environment and obtain expression, metadata and graph inputs.
2. [Run full scKNIFE](tutorials/10000-setup.md): sample cells, prepare expression, tune on a representative subset, fit the graph model, reconstruct activities and save the fitted objects.
3. [Cell-type retention](tutorials/10000-annotation.md): compare H-space clusters and source cell types on the same UMAP.
4. [Metabolic-task analysis](tutorials/10000-metabolic.md): calculate MEAI, inspect one task distribution, and test metabolic-task activity.
5. [GO differential activity](tutorials/10000-go.md): compare GO biological-process activities between conditions within cell types.
6. [Factor interpretation](tutorials/10000-factors.md): group the fitted factors into four consensus clusters (CCs), interpret representative node groups (NGs), and examine their activities.
7. [Cell–cell signaling](tutorials/10000-signaling.md): summarize ligand–receptor scores, permutation FDR, network roles, VEGF signaling and descriptive On-minus-Pre differences.

**Optional extension:** [Rank sensitivity and projection](tutorials/10000-validation.md) examines an 8,000/2,000 training–held-out split within the sample. It is separate from the main biological workflow.

## Run as a notebook

Download the {download}`10,000-cell notebook <downloads/scKNIFE_10000_tutorial.ipynb>` to explore the workflow in Jupyter. To execute it from the command line, use `python scripts/run_tutorial.py`; the completed notebook is saved in your output directory.

## Shared inputs and terminology

Metabolic, GO, factor and signaling analyses use the saved fit and activity scores.

- **Source cell type** means the original `cellType` metadata used to group and describe cells. Cell-type retention concerns biological organization in the learned representation.
- **H** contains cell-factor activities; **W** contains biological-node loadings. The factor page uses **CC** for a consensus cluster of factors and **NG** for its selected biological nodes. The fitted rank is selected by pre-fit tuning; the interpretation resolution is k=4.
- **On minus Pre** means on-treatment minus pre-treatment. Metabolic and GO comparisons use exploratory cell-level tests; signaling contrasts show pooled differences. Neither accounts for patient pairing.

Click a figure to view it at full resolution.

## Extend the activity-based analysis

[Other activity-based analyses](tutorials/10000-states.md) illustrates cell-state programs as an optional task. Choose activities for entities represented in the fitted graph. If a new question requires entities or relationships absent from that graph, supply an appropriate graph and refit before reconstructing their activities. Additional metadata or interaction resources may be needed for the chosen analysis.

