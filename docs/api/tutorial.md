# Tutorial functions

Use these functions with the saved 36,721-cell model. Score and differential-analysis functions return tables for filtering and export.

After `import scknife`, use `scknife.tutorial` for setup, `scknife.fp` for factor interpretation, `scknife.activity` for metabolic and GO comparisons, and `scknife.communication` for signaling. These names expose the same functions documented below.

## Setup and clustering

```{eval-rst}
.. autofunction:: scknife.tutorial.core.output_dir

.. autofunction:: scknife.tutorial.core.marker_graph

.. autofunction:: scknife.tutorial.core.build_graph

.. autofunction:: scknife.tutorial.core.optimize_hyperparameters

.. autofunction:: scknife.tutorial.core.fit_and_score

.. autofunction:: scknife.tutorial.core.load_graph

.. autofunction:: scknife.tutorial.core.load_fit

.. autofunction:: scknife.tutorial.core.load_activities

.. autofunction:: scknife.tutorial.core.save_fit

.. autofunction:: scknife.tutorial.core.figure

.. autofunction:: scknife.tutorial.core.sample_bassez

.. autofunction:: scknife.tutorial.core.plot_sample

.. autofunction:: scknife.tutorial.core.prepare_tuning

.. autofunction:: scknife.tutorial.core.select_parameters

.. autofunction:: scknife.tutorial.core.plot_convergence

.. autofunction:: scknife.tutorial.core.plot_clusters

```

## Factor interpretation

```{eval-rst}
.. autofunction:: scknife.tutorial.factors.group_factors

.. autofunction:: scknife.tutorial.factors.plot_usage

.. autofunction:: scknife.tutorial.factors.plot_group_selection

.. autofunction:: scknife.tutorial.factors.plot_consensus

.. autofunction:: scknife.tutorial.factors.plot_resolutions

.. autofunction:: scknife.tutorial.factors.plot_node_dictionary

.. autofunction:: scknife.tutorial.factors.plot_composition

.. autofunction:: scknife.tutorial.factors.plot_entities

.. autofunction:: scknife.tutorial.factors.plot_group_activity

.. autofunction:: scknife.tutorial.factors.plot_associations

.. autofunction:: scknife.tutorial.factors.plot_cell_profiles

.. autofunction:: scknife.tutorial.factors.plot_coverage

.. autofunction:: scknife.tutorial.factors.plot_constituents

```

## Metabolic tasks and GO

```{eval-rst}
.. autofunction:: scknife.tutorial.activity.node_scores

.. autofunction:: scknife.tutorial.activity.calculate_meai

.. autofunction:: scknife.tutorial.activity.plot_task_umap

.. autofunction:: scknife.tutorial.activity.plot_task_violin

.. autofunction:: scknife.tutorial.activity.compare_activity

.. autofunction:: scknife.tutorial.activity.select_results

.. autofunction:: scknife.tutorial.activity.plot_effects

.. autofunction:: scknife.tutorial.activity.plot_volcano

```

## Signaling

```{eval-rst}
.. autofunction:: scknife.tutorial.signaling.score_signaling

.. autofunction:: scknife.tutorial.signaling.plot_matrix

.. autofunction:: scknife.tutorial.signaling.plot_contributions

.. autofunction:: scknife.tutorial.signaling.permutation_fdr

.. autofunction:: scknife.tutorial.signaling.plot_fdr

.. autofunction:: scknife.tutorial.signaling.centrality

.. autofunction:: scknife.tutorial.signaling.plot_centrality

.. autofunction:: scknife.tutorial.signaling.pathway_scores

.. autofunction:: scknife.tutorial.signaling.plot_strength

.. autofunction:: scknife.tutorial.signaling.plot_network

.. autofunction:: scknife.tutorial.signaling.plot_pathway

.. autofunction:: scknife.tutorial.signaling.plot_pathway_pairs

```

## Optional analyses

```{eval-rst}
.. autofunction:: scknife.tutorial.optional.rank_projection

.. autofunction:: scknife.tutorial.optional.plot_projection

.. autofunction:: scknife.tutorial.optional.state_programs

.. autofunction:: scknife.tutorial.optional.plot_states

```
