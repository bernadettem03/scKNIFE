"""Check biological labels and cell ordering exposed by the tutorial API."""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from matplotlib import pyplot as plt
from scknife.tutorial import plot_task_umap
from scknife.tutorial.activity import plot_effects


def test_umap_uses_requested_scores_on_the_original_cells(tmp_path, monkeypatch):
    monkeypatch.setenv("SCKNIFE_TUTORIAL_OUTPUT", str(tmp_path))
    fit = AnnData(np.zeros((3, 1)))
    fit.obsm["X_umap"] = np.array([[0., 1.], [2., 3.], [4., 5.]])
    scores = pd.DataFrame({"one": [3., 1., 2.], "two": [1., 3., 2.]}, index=fit.obs_names)
    ax = plot_task_umap(fit, scores, task="one")
    np.testing.assert_equal(ax.collections[0].get_array().data, [1., 2., 3.])
    np.testing.assert_equal(ax.collections[0].get_offsets().data, fit.obsm["X_umap"][[1, 2, 0]])
    with pytest.raises(ValueError, match="cell order"):
        plot_task_umap(fit, scores.iloc[::-1], task="one")
    plt.close(ax.figure)


def test_reversed_contrast_changes_plot_labels(tmp_path, monkeypatch):
    monkeypatch.setenv("SCKNIFE_TUTORIAL_OUTPUT", str(tmp_path))
    result = pd.DataFrame({"cell_type": ["T_cell"], "cohens_d": [.6], "passes_cutoffs": [True]})
    result.attrs.update(kind="GO_term", score_name="GO-node activity", expansion="E",
                        contrast=("On", "Pre"), min_d=.3)
    ax = plot_effects(result)
    assert "Pre minus On" in ax.get_xlabel()
    assert [item.get_text() for item in ax.get_legend().texts][:2] == ["Higher in On", "Higher in Pre"]
    plt.close(ax.figure)


def test_effects_override_uses_all_cutoffs_without_changing_results(tmp_path, monkeypatch):
    monkeypatch.setenv("SCKNIFE_TUTORIAL_OUTPUT", str(tmp_path))
    result = pd.DataFrame({"cell_type": ["T_cell"] * 4, "cohens_d": [-.5, .49, .7, .8],
                           "log2FC": [-.25, .4, .24, .5], "adj_p_value": [.009, .001, .001, .01],
                           "passes_cutoffs": [False, True, True, True]})
    result.attrs.update(kind="task", score_name="MEAI", expansion="E", contrast=("Pre", "On"),
                        q_cutoff=.05, min_d=.3, min_log2fc=.25)
    ax = plot_effects(result, q_cutoff=.01, min_d=.5)
    from matplotlib.colors import to_hex
    assert [to_hex(c) for c in ax.collections[0].get_facecolors()] == ["#327e9b"] + ["#c4cbd1"] * 3
    assert [line.get_xdata()[0] for line in ax.lines] == [-.5, 0, .5]
    assert result.passes_cutoffs.tolist() == [False, True, True, True]
    plt.close(ax.figure)
