# Cell–cell signaling with 10,000 cells

Follow the dedicated {doc}`tutorials/10000-signaling` page for cell signaling, with its code, plots and interpretation together.

The primary signaling results use refined annotations. The {doc}`tutorial-10000` includes shared-scale sender–receiver heatmaps and ligand–receptor bubbles, with refined results left and source references right.

Prepare the directional interaction CSV using {doc}`tutorial-inputs`. After fitting and both annotation stages, run:

```bash
python scripts/run_tutorial_10000_signaling.py
python scripts/run_tutorial_10000_refined_downstream.py
```

The first script prepares the reference comparison; the second generates primary refined results and matched comparisons under `tutorial/outputs_10000/annotation_downstream_comparison/`. Read `refined_communication_matrix.csv` and `refined_interaction_details.csv.gz` for primary results.

These are pooled, descriptive activity-based scores, not CellChat probabilities or significance tests. No treatment comparison is performed here. Group-based specificity normalization changes with the set of groups; score differences reflect both membership and this normalization. Uncertain is not a biological type. Inspect rare-group counts and component coverage before interpreting scores. Complete executed code and current numerical checks are in the tutorial.
