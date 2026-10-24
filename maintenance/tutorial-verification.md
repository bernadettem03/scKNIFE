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

