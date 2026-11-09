"""Check that the displayed table and volcano use the same cells and cutoffs."""
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from matplotlib import pyplot as plt
from matplotlib.colors import to_hex
from scknife import activity


def test_custom_metadata_and_display_choices(tmp_path, monkeypatch):
    monkeypatch.setenv("SCKNIFE_TUTORIAL_OUTPUT", str(tmp_path))
    obs = pd.DataFrame({
        "lineage": np.repeat(["A", "B"], 20),
        "visit": np.tile(np.repeat(["Before", "After"], 10), 2),
        "batch": np.tile(["one", "two"], 20),
    }, index=[f"cell{i}" for i in range(40)])
    fit = AnnData(obs=obs)
    change = np.where(obs.visit == "After", 3., 0.)
    baseline = 1 + np.arange(len(obs)) % 5 / 10
    scores = pd.DataFrame({"up": baseline + change, "down": baseline + 3 - change,
                           "steady": baseline}, index=obs.index)
    settings = dict(cell_type_key="lineage", condition_key="visit",
                    contrast=("Before", "After"), min_cells=5)
    result = activity.compare_activity(fit, scores, kind="task", score_name="MEAI",
                                       subset=None, **settings)
    assert (result[["n_group1", "n_group2"]] == 10).all().all()
    saved = result.copy()
    for cell_type in ["A", "B"]:
        for q in [.05, 0.]:
            cutoffs = dict(q_cutoff=q, min_d=.5, min_log2fc=.25)
            passing = activity.select_results(result, cell_type=cell_type, **cutoffs)
            ax = activity.plot_volcano(result, cell_type=cell_type, **cutoffs)
            rows = result.loc[result.cell_type.eq(cell_type)]
            colored = [to_hex(color) != "#c4cbd1" for color in ax.collections[0].get_facecolors()]
            assert set(rows.loc[colored, "feature"]) == set(passing.feature)
            assert set(passing.feature) == ({"up", "down"} if q else set())
            assert passing.adj_p_value.is_monotonic_increasing
            plt.close(ax.figure)
    pd.testing.assert_frame_equal(result, saved)
    with pytest.raises(ValueError, match="not tested"):
        activity.select_results(result, cell_type="unknown")
    filtered = activity.compare_activity(fit, scores, kind="task", score_name="MEAI",
                                         subset={"batch": "one"}, **settings)
    assert (filtered[["n_group1", "n_group2"]] == 5).all().all()
    ax = activity.plot_task_violin(fit, scores, "up", cell_type="A", subset={"batch": "one"}, **settings)
    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["Before", "After"]
    plt.close(ax.figure)
    fit.obs["expansion"] = obs.batch.map({"one": "E", "two": "NE"})
    legacy = activity.compare_activity(fit, scores, kind="task", score_name="MEAI",
                                       expansion="E", **settings)
    pd.testing.assert_frame_equal(legacy, filtered)
