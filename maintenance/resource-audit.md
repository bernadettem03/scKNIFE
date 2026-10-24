# Detailed resource audit

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

