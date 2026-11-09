# Metabolic tasks

Explore task scores on the UMAP and compare metabolic activity between conditions.

## Calculate MEAI

```python
import scknife

fit = scknife.tutorial.load_fit()
graph = scknife.tutorial.load_graph()
activity = scknife.tutorial.load_activities(fit)
scores = scknife.activity.calculate_meai(activity, graph)
score_name = "MEAI"
# scores, score_name = scknife.activity.node_scores(activity, graph, "task"), "Task-node activity"
```

MEAI weights task activity by reaction and gene support, then normalizes across tasks within each cell. The two small [human task tables](https://github.com/earmingol/scCellFie/tree/c66c6be030a8ac871d998f17e592477b97781bf9/task_data/homo_sapiens) download on first use.

(test-differential-metabolic-task-activity)=
## Run the differential activity test

Use the metadata columns for cell type and condition. Here `subset={"expansion": "E"}` selects expanded cells in this dataset. Use `subset=None` to include all cells. With `("Pre", "On")`, positive effects indicate higher activity during treatment. The cell-level Wilcoxon test requires 20 cells per condition and treats cells from the same patient as independent.

```python
comparison = dict(
    cell_type_key="cellType", condition_key="timepoint",
    subset={"expansion": "E"}, contrast=("Pre", "On"), min_cells=20,
)
task_de = scknife.activity.compare_activity(
    fit, scores, kind="task", score_name=score_name, **comparison
)
```

(review-the-differential-results)=
## Inspect the results table

Choose a tested cell type and display cutoffs. `passing` contains every highlighted row for its volcano plot, ordered by Benjamini–Hochberg adjusted p-value. Preview the first five.

```python
cell_type = "Cancer_cell"  # Other options: "T_cell", "Myeloid_cell"
cutoffs = dict(q_cutoff=0.05, min_d=0.3, min_log2fc=0.25)
passing = scknife.activity.select_results(task_de, cell_type=cell_type, **cutoffs)
columns = ["cell_type", "feature", "log2FC", "cohens_d", "adj_p_value"]
passing[columns].head().style.format({"log2FC": "{:.3f}", "cohens_d": "{:.3f}", "adj_p_value": "{:.2e}"})
```

```{csv-table}
:file: ../assets/tutorial-36721/metabolic_walkthrough/top_differential_tasks.csv
:header-rows: 1
```

Across all cell types, 40 of 1,526 comparisons pass the example cutoffs.

## Plot the differential results

```python
scknife.activity.plot_effects(task_de, **cutoffs)
```

```{image} ../assets/tutorial-36721/metabolic_walkthrough/task_effects.png
:class: compact-scientific-plot analysis-single
:align: center
:alt: Positive effects indicate higher On-treatment activity
```

Each point shows one task’s standardized effect within a cell type. Colored points have adjusted p < 0.05, |Cohen’s d| ≥ 0.3 and |log₂ fold change| ≥ 0.25. Edit `cutoffs` to update the highlights without rerunning the test.

### Inspect one cell type

```python
scknife.activity.plot_volcano(task_de, cell_type=cell_type, **cutoffs)
```

```{image} ../assets/tutorial-36721/metabolic_walkthrough/task_volcano.png
:class: compact-scientific-plot analysis-wide
:align: center
:alt: Each point is a tested metabolic task
```

Each point is a tested metabolic task. The guide lines mark the effect-size and adjusted p-value cutoffs. Colored points also pass the absolute log₂ fold-change cutoff. Change `cell_type` or `cutoffs` above to explore the same results.

### View task activity on the UMAP

```python
tasks = ["Glucose to lactate conversion", "GMP salvage from guanine"]
# For one panel: tasks = "Mevalonate synthesis"
scknife.tutorial.plot_task_umap(fit, scores, task=tasks, score_name=score_name)
```

```{image} ../assets/tutorial-36721/metabolic_walkthrough/meai_task_umaps.png
:class: compact-scientific-plot analysis-wide
:align: center
:alt: Each panel shows the same cells and coordinates
```

Each panel shows the same cells and coordinates. Brighter colors indicate higher scores within that task; the color limits are separate.

### Compare one task

```python
scknife.activity.plot_task_violin(fit, scores, task="GMP salvage from guanine",
                     cell_type=cell_type, score_name=score_name, **comparison)
```

```{image} ../assets/tutorial-36721/metabolic_walkthrough/metabolic_task_distribution.png
:class: compact-scientific-plot metabolic-violin
:align: center
:alt: The white dots mark medians and the dark bars the middle 50% of scores
```

The white dots mark medians and the dark bars the middle 50% of scores. The example compares cancer cells in the expanded group.


{download}`Download this example <../downloads/metabolic_walkthrough.py>`.
