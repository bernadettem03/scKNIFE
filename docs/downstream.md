# Downstream analysis of 10,000 cells

Follow the dedicated {doc}`tutorials/10000-metabolic` page for metabolic analysis and paired patient comparisons, with its code, plots and interpretation together.

Primary analysis uses `scknife_breast_refined`, produced by breast-context cluster annotation followed by fixed rare-cell marker gates. See {doc}`tutorial-10000` for executed code, primary plots and matched source-reference comparisons.

Run `scripts/run_tutorial_10000_refined_downstream.py` after fitting, both annotation stages and reference signaling. The runner regenerates state/task means, all 16 H-factor summaries, patient-paired task effects and descriptive signaling. The same 10,000 cells and fitted activities are used for both label sets; graphs, factors and embeddings remain unchanged during regrouping.

Primary refined plots are left and references right, with common features and scales. Patients require ten cells per timepoint and groups need three paired patients; eligibility is recalculated for each annotation. Rare groups and Uncertain remain visible but should be interpreted with their cell counts. Label-independent W node programs are generated from the new fit. See the full tutorial for numerical findings and limitations rather than interpreting reference agreement as independent biological validation.
