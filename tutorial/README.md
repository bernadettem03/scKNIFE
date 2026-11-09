# Tutorial for the 36,721-cell TNBC dataset

[Open the notebook](scKNIFE_Bassez_Cohort_B_complete_tutorial.ipynb) or [start with installation](../docs/installation.md). Place the counts and metadata for the dataset in `tutorial/inputs/` as shown in [Run scKNIFE](../docs/tutorials/36721-setup.md).

```bash
python scripts/run_tutorial.py
```

The notebook uses Version 1, tunes on 2,000 cells, then fits all 36,721 TNBC cells. After clustering, choose factor interpretation, metabolic tasks, GO analysis or cell–cell signaling. Results are saved in `tutorial/outputs_36721/` by default. The synthetic regression fixture is separate.
