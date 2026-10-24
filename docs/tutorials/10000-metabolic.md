# Metabolic and GO differential activity

[Tutorial overview](../tutorial-10000.md) · [Previous: Cell-type retention](10000-annotation.md) · [Next: Factor interpretation](10000-factors.md)

Calculate metabolic scores, inspect one task, then compare metabolic-task and GO activities between conditions.

**Input:** the saved [10,000-cell fit](10000-setup.md). [Download the runnable script](../downloads/metabolic_walkthrough.py).

## Load activities and choose a comparison

Use **On minus Pre** within the expanded (`E`) group. The commented options select non-expanded cells or another cell type. Each cell type needs at least 20 cells in each condition. These are exploratory cell-level comparisons; the tests do not account for cells nested within patients.

```python
from pathlib import Path
import os, json, hashlib, textwrap
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scknife import run_scanpy_KNIFE_DE

BASE = Path(os.environ.get('SCKNIFE_TUTORIAL_OUTPUT', 'tutorial/outputs_10000'))
OUT = BASE / 'metabolic_walkthrough'
OUT.mkdir(parents=True, exist_ok=True)
fit = sc.read_h5ad(BASE / 'fitted_subset.h5ad')
activity = pd.read_parquet(BASE / 'activity_scores.parquet')
assert fit.n_obs == 10000 and fit.obs_names.is_unique
assert activity.columns.equals(fit.obs_names) and activity.index.is_unique
meta = fit.obs[['cellType', 'timepoint', 'expansion']].astype(str)
G = nx.read_graphml('graphs/version_1/unified_graph.graphml')

def nodes_of_kind(kind):
    return sorted(n for n, d in G.nodes(data=True)
                  if d.get('kind') == kind and n in activity.index)

EXPANSION = 'E'
# EXPANSION = 'NE'  # Non-expanded group.
CELL_TYPE = 'Cancer_cell'
# CELL_TYPE = 'T_cell'
# CELL_TYPE = 'Myeloid_cell'
CONTRAST = ('Pre', 'On')  # group 2 minus group 1.
MIN_CELLS = 20
Q_CUTOFF, D_CUTOFF, LOGFC_CUTOFF = 0.05, 0.3, 0.25
GROUP_NAME = {'E': 'Expanded', 'NE': 'Non-expanded'}[EXPANSION]
selected_meta = meta.loc[meta.expansion.eq(EXPANSION)
                         & meta.timepoint.isin(CONTRAST)].copy()
counts = pd.crosstab(selected_meta.cellType, selected_meta.timepoint)
counts = counts.reindex(columns=list(CONTRAST), fill_value=0)
counts['eligible'] = counts[list(CONTRAST)].min(axis=1).ge(MIN_CELLS)
counts.to_csv(OUT / 'condition_cell_counts.csv')
eligible = counts.index[counts.eligible]
test_meta = selected_meta.loc[selected_meta.cellType.isin(eligible)].copy()
assert CELL_TYPE in eligible, 'Choose a cell type with enough cells in both conditions.'
print(counts.to_string())

plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12,
                     'axes.labelsize': 11, 'xtick.labelsize': 10,
                     'ytick.labelsize': 10, 'axes.spines.top': False,
                     'axes.spines.right': False})
def save(name):
    plt.gcf().savefig(OUT / (name + '.png'), dpi=170, bbox_inches='tight')
    plt.close()
```

In this example, 6,924 cells belong to the expanded group. Six cell types, containing 6,875 cells, meet the minimum counts; mast cells and pDCs have too few Pre/On cells for these tests.

## Metabolic-task analysis

### 1. Calculate MEAI

For each task, combine its activity with the mean activity of its linked reactions and genes, then divide by the sum across tasks in that cell:

**MEAI = task activity × √(reaction support × gene support), normalized across tasks per cell.**

Reaction and gene support use the pinned [scCellFie human task tables](https://github.com/earmingol/scCellFie/tree/c66c6be030a8ac871d998f17e592477b97781bf9/task_data/homo_sapiens). The code downloads and verifies these tables once; the support table records the nodes included. MEAI describes relative metabolic support, not metabolic flux. Uncomment the alternative to use raw task-node activity; its scale and differential results can differ.

```python
from urllib.request import urlopen
TASK_TABLES = BASE / 'task_support'
TASK_TABLES.mkdir(exist_ok=True)
TABLE_BASE = ('https://raw.githubusercontent.com/earmingol/scCellFie/'
              'c66c6be030a8ac871d998f17e592477b97781bf9/task_data/homo_sapiens/')
TABLE_HASHES = {
    'Task_by_Rxn.csv': 'dc519e793a83a3fe98874e8ae6636b83d2b12bc3c25b42bd094dd408ba923209',
    'Task_by_Gene.csv': '79059bc0ceec00763e9fa568dc7e2d4fc8f0abc29582da29abd08910a8bf5e66',
}
for filename, expected in TABLE_HASHES.items():
    path = TASK_TABLES / filename
    data = path.read_bytes() if path.exists() else urlopen(TABLE_BASE + filename, timeout=60).read()
    assert hashlib.sha256(data).hexdigest() == expected, f'Unexpected support table: {filename}'
    if not path.exists():
        path.write_bytes(data)
task_rxn = pd.read_csv(TASK_TABLES / 'Task_by_Rxn.csv', index_col=0)
task_gene = pd.read_csv(TASK_TABLES / 'Task_by_Gene.csv', index_col=0)
tasks = nodes_of_kind('task')
raw_task_activity = activity.loc[tasks].T
EPS = 1e-10
support_rows, combined = [], []
for task in tasks:
    reactions = task_rxn.columns[task_rxn.loc[task].eq(1)].tolist()
    genes = task_gene.columns[task_gene.loc[task].eq(1)].tolist()
    assert set(reactions + genes).issubset(activity.index), f'Missing activity nodes for {task}'
    assert reactions and genes, f'Missing scored support for {task}'
    reaction_support = activity.loc[reactions].mean(axis=0).to_numpy()
    gene_support = activity.loc[genes].mean(axis=0).to_numpy()
    task_activity = activity.loc[task].to_numpy()
    combined.append((task_activity + EPS)
                    * np.sqrt((reaction_support + EPS) * (gene_support + EPS)))
    support_rows.append({'task': task, 'reactions': len(reactions), 'genes': len(genes),
                         'reaction_nodes': ';'.join(reactions), 'gene_nodes': ';'.join(genes)})
combined = np.column_stack(combined)
meai = pd.DataFrame(combined / combined.sum(axis=1, keepdims=True),
                    index=fit.obs_names, columns=tasks)
assert np.isfinite(meai.to_numpy()).all() and (meai.to_numpy() >= 0).all()
assert np.allclose(meai.sum(axis=1), 1)
pd.DataFrame(support_rows).to_csv(OUT / 'meai_support.csv', index=False)
meai.to_parquet(OUT / 'meai_scores.parquet')

metabolic_scores, SCORE_NAME = meai, 'MEAI'
# metabolic_scores, SCORE_NAME = raw_task_activity, 'Task-node activity'
```

### 2. Inspect one metabolic task

Compare the cell-level distributions for one task in the selected cell type and expansion group.

```python
EXAMPLE_TASK = 'ATP generation from glucose (hypoxic conditions) - glycolysis'
# EXAMPLE_TASK = 'ATP regeneration from glucose (normoxic conditions) - glycolysis + krebs cycle'
# EXAMPLE_TASK = 'Acetoacetate synthesis'
example_meta = test_meta.loc[test_meta.cellType.eq(CELL_TYPE)]
values = [metabolic_scores.loc[example_meta.index[example_meta.timepoint.eq(c)],
                               EXAMPLE_TASK].to_numpy() for c in CONTRAST]
fig, ax = plt.subplots(figsize=(6.5, 4.2))
violins = ax.violinplot(values, showmedians=True, showextrema=False)
for body, color in zip(violins['bodies'], ['#327e9b', '#c86556']):
    body.set_facecolor(color); body.set_alpha(0.65)
violins['cmedians'].set_color('#273746')
ax.set_xticks([1, 2], [f'{c}\n(n={len(v):,})' for c, v in zip(CONTRAST, values)])
ax.set_ylabel(SCORE_NAME)
ax.set_title(f'{CELL_TYPE} · {GROUP_NAME}\n' + textwrap.fill(EXAMPLE_TASK, 48))
fig.tight_layout(); save('task_violin')
```

```{image} ../assets/tutorial-10000/metabolic_walkthrough/task_violin.png
:class: compact-scientific-plot
:align: center
:alt: Metabolic score distributions before and during treatment
```

The violin width shows the density of cell scores; the horizontal line marks the median. Each condition is labelled with its cell count. The plotted score follows the choice above.

### 3. Test metabolic-task activity

Run two-sided Wilcoxon rank-sum tests within each eligible cell type. The package computes Cohen’s d, log₂(mean On / mean Pre), and BH-adjusted p-values across all tested task–cell-type combinations in this expansion group. Tests use the selected scores without a log transformation.

Highlight results meeting **q < 0.05**, **|d| ≥ 0.3** and **|log₂FC| ≥ 0.25**. Positive effects indicate higher scores in On.

```python
def differential_activity(scores):
    aligned = scores.loc[test_meta.index]
    assert aligned.index.equals(test_meta.index)
    assert np.isfinite(aligned.to_numpy()).all() and (aligned.to_numpy() >= 0).all()
    import warnings
    with warnings.catch_warnings():
        # Compatibility notices from the pinned Scanpy/Pandas versions.
        warnings.filterwarnings('ignore', message='The behavior of DataFrame.sum with axis=None', category=FutureWarning)
        warnings.filterwarnings('ignore', message='The previous implementation of stack is deprecated', category=FutureWarning)
        result = run_scanpy_KNIFE_DE(
            aligned, test_meta, cell_type_key='cellType', condition_key='timepoint',
            contrasts=[CONTRAST], features=aligned.columns.tolist(), min_cells=MIN_CELLS)
    assert not result.empty
    result['passes_cutoffs'] = (result.adj_p_value.lt(Q_CUTOFF)
                                & result.cohens_d.abs().ge(D_CUTOFF)
                                & result.log2FC.abs().ge(LOGFC_CUTOFF))
    return result

task_de = differential_activity(metabolic_scores)
task_de['score'] = SCORE_NAME
task_de.to_csv(OUT / 'metabolic_differential_activity.csv', index=False)

def effect_overview(result, title, filename):
    fig, ax = plt.subplots(figsize=(8, 4.8))
    rng = np.random.default_rng(42)
    types = sorted(result.cell_type.unique())
    for i, celltype in enumerate(types):
        rows = result.loc[result.cell_type.eq(celltype)]
        colors = np.where(rows.passes_cutoffs,
                          np.where(rows.cohens_d.ge(0), '#c86556', '#327e9b'), '#c4cbd1')
        ax.scatter(rows.cohens_d, i + rng.uniform(-0.22, 0.22, len(rows)),
                   c=colors, s=9, alpha=0.6, linewidths=0, rasterized=True)
    ax.set_yticks(range(len(types)), types)
    ax.invert_yaxis()
    for x in [-D_CUTOFF, 0, D_CUTOFF]:
        ax.axvline(x, color='#526171', lw=0.7, ls='--' if x else '-')
    ax.set_xlabel("Cohen's d · On minus Pre")
    ax.set_title(f'{title} · {GROUP_NAME}')
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=c, marker='o', linestyle='', label=t)
                       for c, t in [('#327e9b', 'Higher in Pre'), ('#c86556', 'Higher in On'),
                                    ('#c4cbd1', 'Below display cutoffs')]],
              loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=3, frameon=False)
    fig.tight_layout(); save(filename)

effect_overview(task_de, f'Metabolic tasks · {SCORE_NAME}', 'task_effects')
```

```{image} ../assets/tutorial-10000/metabolic_walkthrough/task_effects.png
:class: compact-scientific-plot
:align: center
:alt: Cell-level metabolic differential activity across cell types
```

Each point is one task tested within a cell type. Colored points meet all three display cutoffs; grey points do not. The horizontal position gives the direction and standardized magnitude of the difference.

### 4. Examine one cell type

Show all tested tasks for the selected cell type. Labels identify up to three tasks meeting the display cutoffs, ranked by |d| × −log₁₀(q).

```python
def volcano(result, title, filename):
    rows = result.loc[result.cell_type.eq(CELL_TYPE)].copy()
    rows['neg_log_q'] = -np.log10(rows.adj_p_value.clip(lower=np.finfo(float).tiny))
    colors = np.where(rows.passes_cutoffs,
                      np.where(rows.cohens_d.ge(0), '#c86556', '#327e9b'), '#c4cbd1')
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.scatter(rows.cohens_d, rows.neg_log_q, c=colors, s=16, alpha=0.65, linewidths=0)
    ax.axhline(-np.log10(Q_CUTOFF), color='#526171', ls='--', lw=0.8)
    for x in [-D_CUTOFF, D_CUTOFF]:
        ax.axvline(x, color='#526171', ls='--', lw=0.8)
    chosen = rows.loc[rows.passes_cutoffs].copy()
    chosen['label_rank'] = chosen.cohens_d.abs() * chosen.neg_log_q
    labels = chosen.sort_values(['label_rank', 'feature'], ascending=[False, True]).head(3)
    for i, row in enumerate(labels.itertuples()):
        ax.annotate(textwrap.fill(row.feature, 30), (row.cohens_d, row.neg_log_q),
                    xytext=(0.03, 0.97 - i * 0.25), textcoords='axes fraction',
                    fontsize=9, va='top', bbox=dict(facecolor='white', alpha=0.85, edgecolor='none'),
                    arrowprops=dict(arrowstyle='-', color='#526171', lw=0.6))
    ax.set_xlabel("Cohen's d · On minus Pre"); ax.set_ylabel('−log₁₀(BH-adjusted p-value)')
    ax.set_title(f'{title}\n{CELL_TYPE} · {GROUP_NAME}')
    ax.margins(x=0.15, y=0.15)
    fig.tight_layout(); save(filename)

volcano(task_de, f'Metabolic tasks · {SCORE_NAME}', 'task_volcano')
```

```{image} ../assets/tutorial-10000/metabolic_walkthrough/task_volcano.png
:class: compact-scientific-plot
:align: center
:alt: Metabolic differential activity volcano plot
```

Each point is a metabolic task. Farther left or right means a larger standardized difference; higher points have smaller adjusted p-values. Colors use the same three cutoffs as the overview, including the fold-change threshold.

## GO biological-process differential activity

### 1. Compare GO activity across cell types

Use reconstructed GO-node activities with the same cells, comparison and thresholds. GO scores are not MEAI. BH correction is applied separately across the tested GO-term–cell-type combinations; this is differential activity, not GO enrichment.

```python
go_scores = activity.loc[nodes_of_kind('GO_term')].T

go_de = differential_activity(go_scores)
go_de['score'] = 'GO-node activity'
go_de.to_csv(OUT / 'go_differential_activity.csv', index=False)
effect_overview(go_de, 'GO biological processes', 'go_effects')
```

```{image} ../assets/tutorial-10000/metabolic_walkthrough/go_effects.png
:class: compact-scientific-plot
:align: center
:alt: Cell-level GO differential activity across cell types
```

Each point is a GO term. Colors identify terms passing the adjusted-p-value, effect-size and fold-change cutoffs within each cell type. GO terms can share genes, so these points do not represent independent biological processes.

### 2. Examine GO differences in one cell type

```python
volcano(go_de, 'GO biological processes', 'go_volcano')
```

```{image} ../assets/tutorial-10000/metabolic_walkthrough/go_volcano.png
:class: compact-scientific-plot
:align: center
:alt: GO differential activity volcano plot
```

The selected cell type uses the same color and axis conventions as the metabolic volcano plot. Labels identify up to three qualifying GO terms; the exported table contains every tested term.

## Save results

The output directory contains MEAI scores, support mappings, condition cell counts, full differential-activity tables and the five figures.

```python
summary = {
    'fit_cells': int(fit.n_obs), 'selected_group': EXPANSION,
    'group_cells': int(len(selected_meta)), 'tested_cells': int(len(test_meta)),
    'eligible_cell_types': eligible.tolist(), 'metabolic_score': SCORE_NAME,
    'tasks': len(tasks), 'go_terms': go_scores.shape[1],
    'task_tests': len(task_de), 'go_tests': len(go_de),
    'task_hits': int(task_de.passes_cutoffs.sum()), 'go_hits': int(go_de.passes_cutoffs.sum()),
    'contrast': list(CONTRAST), 'minimum_cells': MIN_CELLS,
    'cutoffs': {'q': Q_CUTOFF, 'absolute_d': D_CUTOFF, 'absolute_log2FC': LOGFC_CUTOFF},
    'cell_id_sha256': hashlib.sha256('\n'.join(fit.obs_names).encode()).hexdigest(),
}
(OUT / 'differential_summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
```

Next: [Factor interpretation](10000-factors.md).
