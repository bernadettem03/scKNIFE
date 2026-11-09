"""Shared statistical and plotting helpers for the differential-activity tutorials."""
import textwrap
import numpy as np
import matplotlib.pyplot as plt
from scknife import run_scanpy_KNIFE_DE


def passes_cutoffs(result, q_cutoff, min_d, min_log2fc):
    """Apply a strict adjusted p-value and inclusive absolute effect thresholds."""
    if not 0 <= q_cutoff <= 1:
        raise ValueError("q_cutoff must be between 0 and 1.")
    if not all(np.isfinite(x) and x >= 0 for x in (min_d, min_log2fc)):
        raise ValueError("min_d and min_log2fc must be finite and nonnegative.")
    return (result.adj_p_value.lt(q_cutoff)
            & result.cohens_d.abs().ge(min_d)
            & result.log2FC.abs().ge(min_log2fc))


def differential_activity(scores, test_meta, contrast, min_cells,
                          q_cutoff=0.05, d_cutoff=0.3, logfc_cutoff=0.25):
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
            contrasts=[contrast], features=aligned.columns.tolist(), min_cells=min_cells)
    assert not result.empty
    result['passes_cutoffs'] = passes_cutoffs(result, q_cutoff, d_cutoff, logfc_cutoff)
    return result


def effect_overview(result, title, filename, group_name, d_cutoff=0.3, *, contrast=("Pre","On"), return_axes=False):
    fig, ax = plt.subplots(figsize=(5.9, 3.8))
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
    for x in [-d_cutoff, 0, d_cutoff]:
        ax.axvline(x, color='#526171', lw=0.7, ls='--' if x else '-')
    ax.set_xlabel(f"Cohen's d · {contrast[1]} minus {contrast[0]}")
    ax.set_title(f'{title} · {group_name}')
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], color=c, marker='o', linestyle='', label=t)
                       for c, t in [('#327e9b', f'Higher in {contrast[0]}'), ('#c86556', f'Higher in {contrast[1]}'),
                                    ('#c4cbd1', 'Below display cutoffs')]],
              loc='upper center', bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=False, fontsize=8,
              columnspacing=1.0, handletextpad=0.3)
    fig.tight_layout()
    fig.savefig(filename, dpi=170, bbox_inches='tight')
    if return_axes:
        return ax
    plt.close(fig)


def volcano(result, title, filename, cell_type, group_name,
            q_cutoff=0.05, d_cutoff=0.3, *, contrast=("Pre","On"), return_axes=False):
    rows = result.loc[result.cell_type.eq(cell_type)].copy()
    rows['neg_log_q'] = -np.log10(rows.adj_p_value.clip(lower=np.finfo(float).tiny))
    colors = np.where(rows.passes_cutoffs,
                      np.where(rows.cohens_d.ge(0), '#c86556', '#327e9b'), '#c4cbd1')
    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    ax.scatter(rows.cohens_d, rows.neg_log_q, c=colors, s=16, alpha=0.65, linewidths=0)
    if q_cutoff > 0:
        ax.axhline(-np.log10(q_cutoff), color='#526171', ls='--', lw=0.8)
    for x in [-d_cutoff, d_cutoff]:
        ax.axvline(x, color='#526171', ls='--', lw=0.8)
    chosen = rows.loc[rows.passes_cutoffs].copy()
    chosen['label_rank'] = chosen.cohens_d.abs() * chosen.neg_log_q
    labels = chosen.sort_values(['label_rank', 'feature'], ascending=[False, True]).head(3)
    for i, row in enumerate(labels.itertuples()):
        ax.annotate(textwrap.fill(row.feature, 30), (row.cohens_d, row.neg_log_q),
                    xytext=(1.04, 0.98 - i * 0.33), textcoords='axes fraction',
                    fontsize=8.5, va='top', bbox=dict(facecolor='white', alpha=0.85, edgecolor='none'),
                    arrowprops=dict(arrowstyle='-', color='#526171', lw=0.6))
    ax.set_xlabel(f"Cohen's d · {contrast[1]} minus {contrast[0]}"); ax.set_ylabel('−log₁₀(BH-adjusted p-value)')
    ax.set_title(f'{title}\n{cell_type} · {group_name}')
    ax.margins(x=0.15, y=0.15)
    fig.subplots_adjust(left=0.14, right=0.63 if len(labels) else 0.96, top=0.85, bottom=0.16)
    fig.savefig(filename, dpi=170, bbox_inches='tight')
    if return_axes:
        return ax
    plt.close(fig)
