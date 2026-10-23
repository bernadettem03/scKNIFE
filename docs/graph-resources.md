# Biological resources and references

scKNIFE graph collections combine biological relationships recorded in their edge-source attributes. Cite scKNIFE together with the resources used in an analysis, and report the selected graph file and checksum. The references below identify the resources; the applicable resource-specific terms govern reuse of the biological data.

## Gene sets and functional annotation

- **Hallmark / MSigDB:** Liberzon et al. (2015), *The Molecular Signatures Database hallmark gene set collection*. [Publication](https://doi.org/10.1016/j.cels.2015.12.004); [MSigDB terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).
- **Gene Ontology:** Ashburner et al. (2000), *Gene ontology: tool for the unification of biology*. [Publication](https://doi.org/10.1038/75556); [GO citation policy, release guidance and CC BY 4.0 terms](https://geneontology.org/docs/go-citation-policy/).
- **Enrichr:** Kuleshov et al. (2016), *Enrichr: a comprehensive gene set enrichment analysis web server 2016 update*. [Publication](https://doi.org/10.1093/nar/gkw377). The frozen unified graphs record the source identifier `GO_Biological_Process_2021`; the focused GO collection also retains the identifier `G_process`.

## Metabolism

- **Human-GEM:** [Model repository and release-specific citation guidance](https://github.com/SysBioChalmers/Human-GEM); [CC BY 4.0 license](https://github.com/SysBioChalmers/Human-GEM/blob/main/LICENSE.md). Source labels occur as `Human_GEM` or `Human-GEM` across the supplied files.
- **Metabolic task and reaction resources:** [scCellFie repository and references](https://github.com/earmingol/scCellFie). The construction workflow uses task and reaction association tables; the corresponding labels include `task_rxn`, `rxn_gene`, and `metabolic_task`. Retain the table and source-model provenance with an analysis.

## Cell identity and signaling

- **CellChat:** Jin et al. (2021), *Inference and analysis of cell-cell communication using CellChat*. [Publication](https://doi.org/10.1038/s41467-021-21246-9); [upstream GPL v3 license](https://github.com/jinworks/CellChat/blob/main/LICENSE).
- **PanglaoDB:** Franzen, Gan and Bjorkegren (2019), *PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data*. [Publication](https://doi.org/10.1093/database/baz046); [resource and citation guidance](https://panglaodb.se/faq.html).
- **CellMarker 2.0:** Hu et al. (2023), *CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data*. [Publication](https://doi.org/10.1093/nar/gkac947); [database](https://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/index.html).
- **CellTypist encyclopedia:** [CellTypist resource](https://www.celltypist.org/). The marker-construction workflow combines an encyclopedia export with PanglaoDB and manually assembled marker lists. Preserve the export identity and citations for the curated lists alongside the analysis.

## Additional relationships retained in focused collections

- **ARCHS4:** Lachmann et al. (2018), *Massive mining of publicly available RNA-seq data from human and mouse*. [Publication](https://doi.org/10.1038/s41467-018-03751-6); [resource guidance](https://archs4.org/help). The recorded source identifier is `ARCHS4 TFs Coexp Edges`.
- **KEGG:** [KEGG resource](https://www.kegg.jp/); [resource terms](https://www.kegg.jp/kegg/legal.html). The focused metabolic-task collection retains 17 edges carrying the `KEGG 2021 Human Pathway Gene Edges` label, including edges carrying more than one source label.
- **Breast-tissue protein interactions:** The supplied collections retain the source identifier `breast-PPI`. Report this identifier and the relevant graph checksum when using these relationships. The input-processing workflow includes interaction-confidence and tissue-expression filters.

## Tutorial expression data

Bassez et al. (2021), *A single-cell map of intratumoral changes during anti-PD1 treatment of patients with breast cancer*. [Publication](https://doi.org/10.1038/s41591-021-01323-8); [Lambrechts lab public processed-count data](https://lambrechtslab.sites.vib.be/en/single-cell).

The current tutorial selects 10,000 TNBC cells from cohort 2 with seed 42, stratifying by authentic patient, timepoint and cell type. Input preparation is described in {doc}`tutorial-inputs`. E/NE denote expanded/non-expanded clonotypes. Resource-specific terms apply to the expression data and associated metadata.

## Verified source records

The following records were checked against the construction inputs on 21 October 2026. Export names identify the files used; a year in a library name and a model timestamp are recorded as such. All supplied graph files retain their existing checksums.

| Resource | Verified input or version | License or provenance record |
|---|---|---|
| Bassez cohort 2 | `1867-counts_cells_cohort2.rds`, converted to `.h5ad`; metadata `1871-BIOKEY_metaData_cohort2_web.csv` | Public Lambrechts processed-count source; the cohort-2 analysis explicitly loads the converted counts and this metadata. |
| scCellFie metabolic tables | Commit `c66c6be030a8ac871d998f17e592477b97781bf9` (`v0.4.5-1-gc66c6be`) | MIT, copyright 2024 Wellcome Sanger Institute; the license and human task/reaction tables match the commit byte-for-byte. |
| GO via Enrichr | `GO_Biological_Process_2021.GO_BP.Gene.edges.csv` and `GO_Biological_Process_2021.nodes.csv` | Identifies the GO Biological Process 2021 library; see the GO citation and terms above. |
| KEGG via Enrichr | `KEGG 2021 Human Pathway Gene Edges.csv` and `KEGG Human Nodes 2021.csv` | Identifies the retained 2021 export; see resource-specific terms above. |
| ARCHS4 via Enrichr | `ARCHS4 TFs Coexp Edges.csv` and `ARCHS4 TFs Coexp Nodes.csv` | Identifies the local coexpression export. |
| Human-GEM | `Human-GEM.xml`, exact match to official **v1.19.0** | CC BY 4.0 verified at that release; Git blob `919a408db8aad28103da95530c4484c9c50b3ff2`. |
| MSigDB/Hallmark | `real_msigdb_processed.csv` | Contains gene-set names, descriptions, source URLs and gene symbols; an originating release identifier remains to be recorded. |
| CellChat | `interactions.csv`, `complex.csv`, `cofactors.csv`, `geneInfo.csv` | Content matches snapshot `592ba49cf4b72130e27eb47a217faf5f84a54f07`; export procedure and retained GPLv3 notice are recorded below. |
| Marker and tissue-interaction inputs | `cellMarker2.0.xlsx`, `cellMarker2.0_cancer_breastTissue.csv`, `encyclopedia_table_celltypist.xlsx`, `breast-PPI_new.txt` | PanglaoDB is retrieved through decoupler with `organism="human", license="academic"`; original export records remain to be completed. |

### Metabolic-table provenance and notice

The human `Task_by_Gene.csv`, `Task_by_Rxn.csv` and `Rxn_by_Gene.csv` files match the scCellFie commit above. The separately stored `Rxn_by_Gene.csv` also matches that checkout's human reaction–gene table exactly.

The generation notebook identifies the [secCellFie source repository](https://github.com/LewisLabUCSD/secCellFie) and Recon 2.2 inputs, including `MT_recon_2_2_entrez_sec.mat` and `taskStructure_met_sec.mat`. This establishes the recorded model lineage; upstream input revisions and notices are tracked separately from the scCellFie software license.

The complete verified scCellFie notice is retained here: {download}`scCellFie MIT license <licenses/scCellFie-MIT.txt>`. This notice covers the identified scCellFie material; resource-specific terms continue to apply to upstream biological inputs.

See {doc}`release-checklist` for the remaining records. Verified filenames, the scCellFie revision and its immediate table provenance are already recorded and do not need to be supplied again.

## Official-source review

Official sources were checked on 21 October 2026. These findings distinguish published terms from the identity of the particular export used in scKNIFE. Current website versions are not assigned retrospectively to older local files.

| Resource | Official-source finding | Remaining evidence |
|---|---|---|
| Human-GEM | The local XML matches the [v1.19.0 file](https://github.com/SysBioChalmers/Human-GEM/blob/v1.19.0/model/Human-GEM.xml) by Git blob hash. The [release license](https://github.com/SysBioChalmers/Human-GEM/blob/v1.19.0/LICENSE.md) is CC BY 4.0. | No further version/license input is needed from the maintainer for this identified model. |
| GO | The [official policy](https://geneontology.org/docs/go-citation-policy/) permits reuse under CC BY 4.0 with attribution and requests release identification. | The precise Enrichr export record remains a provenance detail; a separate GO permission letter is not the identified gap. |
| MSigDB/Hallmark | [Official terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp) specify CC BY 4.0 for v6.0–v7.5.1 and v2022.1 onward, with exceptions for designated collections. Hallmark is not among the listed exceptions. Earlier releases have different terms. | Identify the release behind the processed local CSV; do not substitute the current website release. |
| CellChat | The [official repository](https://github.com/jinworks/CellChat/blob/main/LICENSE) supplies GPLv3. | The matching snapshot and export procedure are recorded below; its GPLv3 notice is retained in this repository. The historical installed package version is not inferred from that match. |
| Bassez | The [lab page](https://lambrechtslab.sites.vib.be/en/single-cell) provides processed counts separately from controlled-access EGA raw reads. | An explicit basis for rehosting the derived subset was not located on the inspected page. Raw-read access conditions were not treated as the processed-count license. |
| KEGG | [Official terms](https://www.kegg.jp/kegg/legal.html) describe academic website use and service-provider licensing. | These terms alone do not establish public redistribution coverage for the retained export. |
| ARCHS4 | The [official help page](https://archs4.org/help) labels expression files CC BY 4.0, while also giving noncommercial-use wording and a commercial-contact instruction. | Clarify applicable terms for the particular TF-coexpression export; expression-file terms were not assumed to cover every derivative. |
| PanglaoDB / CellMarker | [PanglaoDB](https://panglaodb.se/faq.html) provides marker-download and citation guidance; [CellMarker 2.0](https://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/index.html) identifies its database and citation. | A clear export-specific redistribution notice was not located on the inspected pages. |
| CellTypist encyclopedia | The [official encyclopedia](https://www.celltypist.org/encyclopedia/Immune/v2) and [source repository](https://github.com/Teichlab/celltypist_wiki) were inspected. The site links to [Sanger terms](https://www.sanger.ac.uk/terms-of-use/), whose general content default is CC BY-NC-ND 4.0 unless otherwise stated. | Establish the terms for the specific spreadsheet export; do not substitute the CellTypist software MIT license or assume general website terms settle the database export. |
| secCellFie | A [GPLv3 notice](https://github.com/LewisLabUCSD/secCellFie/blob/master/CellFie/LICENSE) exists in the upstream `CellFie` subfolder. | Its scope for the exact model/task inputs and their revision remains to be documented; it is not a root-wide license determination. |
| Breast-PPI | The local filename does not establish an upstream publisher. | Original source URL/publication and terms are still needed; no source has been guessed. |

For Human-GEM, the matching model occurs at release commit `022ed5c2c971e27f739e149c4bbe1328de38337f`. This is a content match, not an inference from the XML timestamp. The model license is also available from the official release link above.

## Database provenance: verified 22 October 2026

**CellChat.** The verified R export procedure exports `CellChatDB.human` to the four CSV paths used by the graph workflow. The interaction, complex and gene-information tables match the official database at [commit 592ba49cf4b72130e27eb47a217faf5f84a54f07](https://github.com/jinworks/CellChat/blob/592ba49cf4b72130e27eb47a217faf5f84a54f07/data/CellChatDB.human.rda) after the same R `write.csv`/`read.csv` round trip, including row and column names. The interaction table has 3,233 rows and 28 columns: 1,920 rows labelled `CellChatDB v1` and 1,313 labelled `CellChatDB v2`. This identifies matching database content, not the historical installed package version or date of download. The historical R session information is unavailable.

The export notebook requests `$cofactors`, whereas the matched database field is `$cofactor`; the existing `cofactors.csv` is empty. This file was not changed. The current focused signaling workflow uses the directional interaction CSV, as documented in the main tutorial. The matched snapshot has the [upstream GPLv3 notice](https://github.com/jinworks/CellChat/blob/592ba49cf4b72130e27eb47a217faf5f84a54f07/LICENSE).

**MSigDB.** Original graph-building code loads `MsigDb_data/real_msigdb_processed.csv`; the adjacent demo notebook only reads that file. A later cell-type review contains `h.all.v2025.1.Hs.symbols.gmt`. Comparing gene memberships gives 49 identical Hallmark sets out of 50; `HALLMARK_G2M_CHECKPOINT` differs. Thus the later version record does not establish the release of the graph input. The original processed-CSV release remains unassigned.

**Breast-PPI.** The located notebooks and package loader consume `breast-PPI_new.txt` and apply interaction-confidence and expression filters. The file begins with interaction/gene-count messages, without an upstream citation or release identifier. No producer or original download source was identified in the inspected scripts and notebooks, so no upstream database name has been assigned.

## Retained upstream notices

The following notices are copied verbatim from the identified upstream sources:

- {download}`CellChat GPLv3 notice <licenses/CellChat-GPL-3.0.txt>` from snapshot `592ba49cf4b72130e27eb47a217faf5f84a54f07`.
- {download}`Human-GEM v1.19.0 license <licenses/Human-GEM-v1.19.0-LICENSE.md>` from the matching model release.
- {download}`scCellFie MIT notice <licenses/scCellFie-MIT.txt>` for the matching task/reaction table checkout.

CellChat database provenance: the human database at the pinned snapshot supplies the interaction definitions exported by `scripts/prepare_tutorial_interactions.R`. The tutorial downloads the upstream database and exports the CSV locally; it does not modify the supplied graphs. The expected interaction-export SHA-256 is `098d8c21cf43d486a548dfc6e4ecdca89221254d6b378ab331973efe98ad4905`. Cite Jin et al. (2021), linked above, with the snapshot identifier. Resource-specific notices accompany the identified upstream material and do not replace the remaining provenance records in the release checklist.
