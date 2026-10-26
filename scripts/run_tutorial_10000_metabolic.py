from pathlib import Path
import os, json, hashlib, textwrap
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scknife.tutorial_differential import differential_activity, effect_overview, volcano

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

plt.rcParams.update({'font.size': 11, 'axes.titlesize': 12,
                     'axes.labelsize': 11, 'xtick.labelsize': 10,
                     'ytick.labelsize': 10, 'axes.spines.top': False,
                     'axes.spines.right': False})
def save(name):
    plt.gcf().savefig(OUT / (name + '.png'), dpi=170, bbox_inches='tight')
    plt.close()


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


UMAP_TASKS = ['Glucose to lactate conversion', 'GMP salvage from guanine']
# UMAP_TASKS = ['Mevalonate synthesis', 'Cholesterol synthesis']
assert 'X_umap' in fit.obsm, 'Use the saved UMAP from the fitted tutorial model.'
xy = np.asarray(fit.obsm['X_umap'])
assert xy.shape == (10000, 2) and np.isfinite(xy).all()
assert metabolic_scores.index.equals(fit.obs_names)
fig, axes = plt.subplots(1, len(UMAP_TASKS), figsize=(3.8 * len(UMAP_TASKS), 3.6),
                         squeeze=False, constrained_layout=True)
for ax, task in zip(axes.ravel(), UMAP_TASKS):
    scores = metabolic_scores[task].to_numpy()
    order = np.argsort(scores, kind='stable')
    points = ax.scatter(xy[order, 0], xy[order, 1], c=scores[order],
                        s=4, cmap='viridis', vmin=0, vmax=scores.max(),
                        linewidths=0, rasterized=True)
    ax.set_title(textwrap.fill(task, 28), fontsize=11)
    ax.set_aspect('equal'); ax.set_axis_off()
    colorbar = fig.colorbar(points, ax=ax, fraction=0.035, pad=0.025, shrink=0.75)
    colorbar.set_label(SCORE_NAME, fontsize=10)
    colorbar.ax.tick_params(labelsize=9)
pd.DataFrame(xy, index=fit.obs_names, columns=['UMAP1', 'UMAP2']).join(
    metabolic_scores[UMAP_TASKS]).to_csv(OUT / 'displayed_task_umaps.csv')
save('meai_task_umaps')


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

EXAMPLE_TASK = 'GMP salvage from guanine'
# EXAMPLE_TASK = 'Glucose to lactate conversion'
# EXAMPLE_TASK = 'Mevalonate synthesis'
example_meta = test_meta.loc[test_meta.cellType.eq(CELL_TYPE)]
values = [metabolic_scores.loc[example_meta.index[example_meta.timepoint.eq(c)],
                               EXAMPLE_TASK].to_numpy() for c in CONTRAST]
positions = [0.8, 1.35]
fig, ax = plt.subplots(figsize=(3.6, 4.2))
violins = ax.violinplot(values, positions=positions, widths=0.28,
                       showmedians=False, showextrema=False)
for body, color in zip(violins['bodies'], ['#327e9b', '#c86556']):
    body.set_facecolor(color); body.set_edgecolor(color)
    body.set_linewidth(0.8); body.set_alpha(0.72)
for position, value in zip(positions, values):
    q1, median, q3 = np.quantile(value, [0.25, 0.5, 0.75])
    ax.vlines(position, q1, q3, color='#243344', linewidth=2.2)
    ax.scatter(position, median, s=20, c='white', edgecolors='#243344',
               linewidths=0.7, zorder=3)
ax.set_xticks(positions, list(CONTRAST))
ax.set_xlim(0.45, 1.7)
ax.set_ylabel(SCORE_NAME)
ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useMathText=True)
ax.set_title(textwrap.fill(EXAMPLE_TASK, 26), fontsize=12, pad=24)
ax.text(0.5, 1.02, f'{CELL_TYPE.replace("_", " ")} · {GROUP_NAME}',
        transform=ax.transAxes, ha='center', fontsize=10, color='#526171')
fig.tight_layout(); save('metabolic_task_distribution')


task_de = differential_activity(metabolic_scores, test_meta, CONTRAST, MIN_CELLS,
                                Q_CUTOFF, D_CUTOFF, LOGFC_CUTOFF)
task_de['score'] = SCORE_NAME
task_de.to_csv(OUT / 'metabolic_differential_activity.csv', index=False)


effect_overview(task_de, f'Metabolic tasks · {SCORE_NAME}',
                OUT / 'task_effects.png', GROUP_NAME, D_CUTOFF)


volcano(task_de, f'Metabolic tasks · {SCORE_NAME}', OUT / 'task_volcano.png',
        CELL_TYPE, GROUP_NAME, Q_CUTOFF, D_CUTOFF)


summary = {
    'fit_cells': int(fit.n_obs), 'selected_group': EXPANSION,
    'group_cells': int(len(selected_meta)), 'tested_cells': int(len(test_meta)),
    'eligible_cell_types': eligible.tolist(), 'metabolic_score': SCORE_NAME,
    'tasks': len(tasks), 'displayed_tasks': UMAP_TASKS, 'violin_task': EXAMPLE_TASK,
    'task_tests': len(task_de),
    'task_hits': int(task_de.passes_cutoffs.sum()),
    'contrast': list(CONTRAST), 'minimum_cells': MIN_CELLS,
    'cutoffs': {'q': Q_CUTOFF, 'absolute_d': D_CUTOFF, 'absolute_log2FC': LOGFC_CUTOFF},
    'cell_id_sha256': hashlib.sha256('\n'.join(fit.obs_names).encode()).hexdigest(),
}
(OUT / 'differential_summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
