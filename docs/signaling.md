# Cell–cell signaling with 10,000 cells

Follow {doc}`tutorials/10000-signaling` for the seven-step signaling analysis. It uses the same fitted 10,000 cells, reconstructed node activities and source `cellType` labels as the main workflow.

Prepare the directional interaction CSV using {doc}`tutorial-inputs`, complete {doc}`tutorials/10000-setup`, then run:

```bash
python scripts/run_tutorial_10000_signaling.py
```

Examples cover sender–receiver activity, leading ligand–receptor pairs, VEGF signaling, permutation FDR, network centrality and descriptive On-minus-Pre contrasts. Outputs are saved under `signaling/` in `SCKNIFE_TUTORIAL_OUTPUT`. The detailed page explains the assumptions and interpretation of each plot. Signaling scores are model-derived summaries, not direct measurements of communication.
