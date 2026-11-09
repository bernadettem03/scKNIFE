# Factor interpretation

The model rank describes how many factors were fitted. Here we group the 24 fitted factors into **four programs** using similarities in their node loadings. The four groups contain 4, 8, 7 and 5 factors, respectively.

```python
import scknife

fit = scknife.tutorial.load_fit()
programs = scknife.fp.group_factors(fit, scknife.tutorial.load_graph(), n_groups=4, bootstraps=200, seed=42)
programs.selected.head()
```

## Factor usage

```python
scknife.fp.plot_usage(programs)
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/factor_usage.png
:class: compact-scientific-plot analysis-wide
:align: center
:alt: Each row is a fitted factor
```

Each row is a fitted factor. Color shows its mean H value within each source cell type.

## Choose a grouping resolution

```python
scknife.fp.plot_group_selection(programs)
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/choose_k.png
:class: compact-scientific-plot analysis-wide
:align: center
:alt: The curves compare ambiguity and separation as the number of programs changes
```

The curves compare ambiguity and separation as the number of programs changes. The dashed line marks the four-group resolution used here.

## Consensus membership

```python
scknife.fp.plot_consensus(programs)
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/consensus.png
:class: compact-scientific-plot analysis-wide
:align: center
:alt: Each tile is the fraction of bootstrap samples in which two factors clustered together
```

Each tile is the fraction of bootstrap samples in which two factors clustered together.

## Compare nearby resolutions

```python
scknife.fp.plot_resolutions(programs)
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/resolution_sensitivity.png
:class: compact-scientific-plot analysis-wide
:align: center
:alt: Compare factor memberships at two, three and four groups
```

Compare factor memberships at two, three and four groups. Group numbers are specific to each row.

## Representative nodes across factors

```python
scknife.fp.plot_node_dictionary(programs, by="factor")
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/node_factors_all.png
:class: compact-scientific-plot combined-entity-plot
:align: center
:alt: Each row is one selected biological node; columns are individual factors ordered by program
```

Each row is one selected biological node; columns are individual factors ordered by program. Compare relative loadings within a row.

## Representative nodes across programs

```python
scknife.fp.plot_node_dictionary(programs, by="group")
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/node_dictionary_all.png
:class: compact-scientific-plot combined-entity-plot
:align: center
:alt: The same nodes are shown against program-mean loadings
```

The same nodes are shown against program-mean loadings. Horizontal lines separate node groups (NGs).

## Node categories

```python
scknife.fp.plot_composition(programs)
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/composition.png
:class: compact-scientific-plot analysis-single
:align: center
:alt: The stacked bars count the kinds of biological nodes represented in each program
```

The stacked bars count the kinds of biological nodes represented in each program.

## One entity category

List the available node categories, factor groups (CC) and corresponding node groups (NG).

```python
programs.entity_kinds, programs.groups, programs.ngroups
```

Display gene sets across all four programs.

```python
# Other kinds: "gene", "GO_term", "reaction", "task", "cellchat_entity"
scknife.fp.plot_entities(programs, kind="geneset")
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/category_all_geneset.png
:class: compact-scientific-plot combined-entity-plot
:align: center
:alt: Each row is a gene set
```

Each row is a gene set. Compare its loading profile across the four programs.

## Node-group activity

```python
scknife.fp.plot_group_activity(programs, group="NG")
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/node_activity.png
:class: compact-scientific-plot analysis-single
:align: center
:alt: Reconstructed node-group activities are summarized by cell type
```

Reconstructed node-group activities are summarized by cell type. Color is standardized within each row.

## Factor-group activity

```python
scknife.fp.plot_group_activity(programs, group="CC")
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/factor_activity.png
:class: compact-scientific-plot analysis-single
:align: center
:alt: Each row summarizes H activity for the factors in a consensus cluster (CC)
```

Each row summarizes H activity for the factors in a consensus cluster (CC).

## Connect nodes, factors and cell types

```python
scknife.fp.plot_associations(programs)
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/associations.png
:class: compact-scientific-plot analysis-wide
:align: center
:alt: Lines connect node groups to factor groups and then to cell types
```

Lines connect node groups to factor groups and then to cell types. Their colors show the sign of the relative association.

## Programs across cells

```python
scknife.fp.plot_cell_profiles(programs)
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/all_cells.png
:class: compact-scientific-plot analysis-wide
:align: center
:alt: Columns are cells ordered by cell type, patient and timepoint
```

Columns are cells ordered by cell type, patient and timepoint. The top tracks identify those groupings.

## Patient representation

```python
scknife.fp.plot_coverage(programs)
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/coverage.png
:class: compact-scientific-plot analysis-single
:align: center
:alt: The heatmap counts cells from each patient within each cell type
```

The heatmap counts cells from each patient within each cell type.

## Selected constituents

```python
scknife.fp.plot_constituents(programs)
```

```{image} ../assets/tutorial-36721/factor_interpretation_k4/constituents.png
:class: compact-scientific-plot analysis-wide
:align: center
:alt: The violins show reconstructed activity for three representative NG1 nodes across cell types
```

The violins show reconstructed activity for three representative NG1 nodes across cell types.

{download}`Download this example <../downloads/factor_interpretation.py>`.
