# Tutorial verification records

## Clean-run verification

The updated command sequence was executed in an isolated output directory using the prepared source inputs. All six stages completed, including 26 tuning fits and 1,000 fresh signaling permutations. The same 10,000 cell identifiers were retained. All 34 main tutorial figures and the optional state-program figure were regenerated; all graph-file checksums remained unchanged.

| Stage | Wall time (seconds) |
| --- | ---: |
| Sampling, tuning, full fit and activity reconstruction | 470.30 |
| Cell-type retention | 2.10 |
| Metabolic analysis | 9.91 |
| Factor interpretation | 12.96 |
| Signaling, including 1,000 permutations | 45.08 |
| Optional rank/projection diagnostic | 117.08 |

Times include process startup and input loading on the tested machine; installation and input preparation are excluded. Tuning accounts for 327 seconds of the first stage. Peak process memory was approximately 3.57 GiB for tuning/full fitting and 2.31 GiB for the optional diagnostic. Other stages were not separately memory-profiled.

Download the {download}`clean-run verification <../docs/downloads/tutorial_clean_run_verification.json>`, {download}`fit report <../docs/downloads/tutorial_10000_run_report.json>` and {download}`projection report <../docs/downloads/tutorial_10000_validation_report.json>`.

CI checks that the six command scripts match the documented Python code, alongside package tests, supplied-graph validation and a strict documentation build. The external-data biological run is verified separately; CI's bundled synthetic 1,000-cell fixture remains a regression check.


## Metabolic and GO differential activity — 24 October 2026

- Executed all 13 code/figure cells of the replacement metabolic section on the existing 10,000-cell fit; five figures generated. Other notebook code and outputs preserved.
- Cell identity hash: `a47ad5d98b131c0e47129c7657dcc35515c675cddd510c60e89f986c21212d63`.
- Expanded group: 6,924 cells; 6,875 cells in six types meet the 20-cell minimum per condition.
- MEAI uses the pinned scCellFie task-support tables, byte-identical to the original analysis inputs. Comparison with the original MEAI function: maximum absolute difference 1.19e-7 (float32 summation order).
- Tested 218 tasks (1,308 comparisons) and 6,036 GO terms (36,216 comparisons). Display cutoffs: q < 0.05, |d| >= 0.3, |log2FC| >= 0.25. Twenty task–cell-type and 405 GO-term–cell-type comparisons pass all three cutoffs.
- Independently verified 12 sampled Wilcoxon p-values, effect sizes and fold changes; recomputed BH corrections across each full result table. The raw task-activity option also completed all 1,308 comparisons.
- Tests are pooled-cell exploratory comparisons, without patient adjustment. No graph or fitting changes.

## Separate GO tutorial page — 24 October 2026

- Added GO as its own page, executable script and notebook stage after metabolic tasks.
- Executed GO first in a fresh Python kernel (6 code/figure cells, 78.8 seconds), then metabolic analysis in another fresh kernel (9 cells, 15.3 seconds). Both use the saved 10,000-cell fit.
- All five figure PNGs, both differential-result CSVs and the MEAI score matrix are byte-identical to the combined-page results.
- Shared statistical/plotting helpers accept explicit settings; GO does not depend on metabolic variables or outputs. The notebook preserves all other stage outputs.
