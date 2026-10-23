# Remaining release records

This checklist contains only items still open after the source audit and 23 October 2026 documentation update. Human-GEM v1.19.0 and its CC BY 4.0 license are now verified by an exact upstream file match and removed from this checklist. Completed source identities and notices are recorded in {doc}`graph-resources`. Existing download notes, license files, source URLs with release information, or permission correspondence can complete these records. A suitable existing license may provide the needed coverage without separate correspondence.

## Resource records to supply or locate

For redistribution coverage, provide an applicable license/terms link, retained download notice, or permission correspondence. For source identification, a release name, original download URL or export record is sufficient to continue the audit. Human-GEM, the matching scCellFie tables and CellChat snapshot do not need to be identified again.

- [ ] **KEGG:** download/export record and applicable redistribution terms for the retained KEGG 2021 relationships. The export year and filenames are already identified.
- [ ] **ARCHS4:** original coexpression-export source/date or version and accompanying terms. The local edge and node filenames are already identified.
- [ ] **Breast-PPI:** original upstream dataset/publication or source URL and applicable terms for `breast-PPI_new.txt`.
- [ ] **Marker resources:** PanglaoDB, CellMarker 2.0 and CellTypist export/version and terms records; sources or authorship notes for the manually curated canonical/TNBC marker lists.
- [ ] **MSigDB/Hallmark:** originating release/download record so the applicable published terms can be attached. Modern-release terms have been found publicly.
- [ ] **GO/Enrichr:** original export/download record for the identified GO Biological Process 2021 library. The official CC BY 4.0 terms and attribution requirements are located; no separate GO permission letter is being requested.
- [ ] **Upstream metabolic inputs:** revision/export record and license scope for the secCellFie/Recon 2.2 inputs referenced by the generation notebook. A GPLv3 notice was located in upstream secCellFie/CellFie; its coverage of these inputs remains to be established. The scCellFie commit, MIT notice and matching task/reaction tables are already verified.

These entries identify records to complete; they do not indicate that permission has been refused. All supplied graphs remain unchanged.

The MSigDB record must identify the older `real_msigdb_processed.csv`; the later 2025.1.Hs Hallmark file differs in one gene set and does not resolve it. For breast-PPI, an original source URL, generating script or the collaborator who supplied it would enable the next trace.

## Completed release preparation

- [x] **Synthetic regression data:** the bundled Bassez H5AD and metadata have been replaced with independently simulated counts and labels. The main 10,000-cell tutorial obtains real inputs from the original source. The public-release candidate begins with a fresh synthetic-only history. Earlier fixture content is preserved solely in separate private archives and backups.

- [x] **Clean release installed:** the new public `bernadettem03/scKNIFE` repository contains the reviewed 239-file release with a fresh synthetic-only root commit, including the updated README and graph guide. Manuscript figure-reproduction folders and the original history are excluded.
- [x] **Original repository preserved:** `bernadettem03/scKNIFE-private-archive-20261023` is private and archived read-only. Separate repository and working-folder backups are retained.
- [x] **Previous commit references checked:** the checked previous commit IDs are unavailable through the new repository; an old commit URL also returned 404 in the browser. The history remains available in the private archive.
- [x] **Working clone synchronized:** the designated server clone matches the new private repository. All 19 supplied graph files are byte-for-byte unchanged.
- [x] **Prior GitHub validation passed:** Documentation and Tests and executable tutorials both completed successfully before fixture replacement. The synthetic working copy passed local regression checks; fresh GitHub checks run on the replacement commit. CI executes the bundled regression workflow; the complete 10,000-cell workflow was tested separately as described below.

## Remaining publication decisions

- [ ] **Complete applicable resource records:** resolve redistribution coverage for bundled third-party inputs and attach applicable notices. Version and export records establish which terms apply; they are not all requests for separate permission letters.
- [x] **Public visibility:** the maintainer authorized publication on 23 October 2026. The clean repository is public; archived repositories remain private. Publication does not resolve the resource provenance records listed above.

## Read the Docs setup

- [x] **Project created:** `scknife-docs`, using the public GitHub repository and branch `main` with `.readthedocs.yaml`.
- [x] **Automatic updates:** a signed GitHub push webhook connects this repository to Read the Docs; SSL verification is enabled. No GitHub App installation is required for this manual integration.
- [x] **Hosted documentation:** the initial documentation is available at [scKNIFE documentation](https://scknife-docs.readthedocs.io/en/latest/).
- [x] **Documentation links:** the hosted URL is included in the README and package metadata.

## Paper citation

- [ ] **Paper citation:** provide the preferred manuscript title, author list and DOI or preprint URL when available. The software citation remains provisional until those details are supplied.

## Verification scope

The 10,000-cell notebook and six analysis stages previously completed in a fresh environment. The candidate preserves their analysis code and outputs; validation results are recorded in the separate review report. Documentation hosting displays the saved notebook outputs; it does not rerun the full biological analysis during a documentation build.

The new repository starts with a fresh root commit; the previous history is preserved only in private backups and the private archive. Subsequent maintenance commits belong to the clean release history.

Synthetic replacement validation: all 19 package tests, the complete regression notebook, regression plotting script, tutorial code-contract checks and strict Sphinx build passed on 23 October 2026. Regeneration reproduced identical counts and metadata. The plotting run completed in 28.5 seconds. These checks were run on the updated working copy; they are separate from the successful GitHub CI on the preceding commit.
