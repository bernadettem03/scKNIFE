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
OUT = BASE / 'go_differential'
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


go_scores = activity.loc[nodes_of_kind('GO_term')].T
go_de = differential_activity(go_scores, test_meta, CONTRAST, MIN_CELLS,
                              Q_CUTOFF, D_CUTOFF, LOGFC_CUTOFF)
go_de['score'] = 'GO-node activity'
go_de.to_csv(OUT / 'go_differential_activity.csv', index=False)
effect_overview(go_de, 'GO biological processes', OUT / 'go_effects.png',
                GROUP_NAME, D_CUTOFF)


volcano(go_de, 'GO biological processes', OUT / 'go_volcano.png',
        CELL_TYPE, GROUP_NAME, Q_CUTOFF, D_CUTOFF)


summary = {
    'fit_cells': int(fit.n_obs), 'selected_group': EXPANSION,
    'group_cells': int(len(selected_meta)), 'tested_cells': int(len(test_meta)),
    'eligible_cell_types': eligible.tolist(), 'score': 'GO-node activity',
    'go_terms': go_scores.shape[1], 'go_tests': len(go_de),
    'go_hits': int(go_de.passes_cutoffs.sum()),
    'contrast': list(CONTRAST), 'minimum_cells': MIN_CELLS,
    'cutoffs': {'q': Q_CUTOFF, 'absolute_d': D_CUTOFF, 'absolute_log2FC': LOGFC_CUTOFF},
    'cell_id_sha256': hashlib.sha256('\n'.join(fit.obs_names).encode()).hexdigest(),
}
(OUT / 'differential_summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
