# Obtain and prepare the 10,000-cell tutorial inputs

The quick start and full tutorial use the same files, sample seed (42), Version 1 graph and output directory. You download the source cohort once; the script selects exactly 10,000 TNBC cells with donor, timepoint and reference cell-type coverage. You do not need to download a separate preselected 10,000-cell dataset.

## 1. Get the Bassez cohort-2 data

Visit the [Lambrechts lab Bassez data page](https://lambrechtslab.sites.vib.be/en/single-cell), scroll to **Data files**, and obtain the **cohort 2** processed counts and matching metadata:

| File | Purpose |
|---|---|
| `1867-counts_cells_cohort2.rds` | Original gene-by-cell count matrix |
| `1871-BIOKEY_metaData_cohort2_web.csv` | Matching cell barcodes, patient IDs, timepoints and cell types |

Choose cohort 2 for both files. If your browser changes a filename, use its actual path in the commands below. The processed counts are the inputs used here; raw EGA sequencing reads are not required. If the lab's file links are unavailable, contact the source lab for these exact processed files. Source terms and citations are recorded in {doc}`graph-resources`.

The original matrix contains 50,693 cells, including 36,721 TNBC cells. Those numbers describe the source pool only. Every main analysis figure is generated after selection of 10,000 cells.

## 2. Convert the R counts to AnnData

Install R with the Matrix package, and complete {doc}`installation` for the Python environment. Convert the source files with:

```bash
mkdir -p tutorial/inputs/converted_counts
Rscript scripts/prepare_bassez_counts.R /absolute/path/1867-counts_cells_cohort2.rds tutorial/inputs/converted_counts
python scripts/prepare_bassez_counts.py tutorial/inputs/converted_counts tutorial/inputs/bassez_cohort2_counts.h5ad
```

The converter transposes genes × cells to cells × genes and preserves values, gene symbols and cell barcodes. Temporary Matrix Market files need extra disk space. If you already have the equivalent unnormalized H5AD file, skip conversion. Do not supply normalized or log-transformed expression: the tutorial performs that step itself.

## 3. Obtain the directional interaction table for signaling

Use the human CellChat database snapshot at commit `592ba49cf4b72130e27eb47a217faf5f84a54f07`. The export below matches the interaction table used for the displayed signaling results. It only needs base R; installing the CellChat package is unnecessary.

```bash
curl --fail --location https://raw.githubusercontent.com/jinworks/CellChat/592ba49cf4b72130e27eb47a217faf5f84a54f07/data/CellChatDB.human.rda --output tutorial/inputs/CellChatDB.human.rda
Rscript scripts/prepare_tutorial_interactions.R tutorial/inputs/CellChatDB.human.rda tutorial/inputs/interactions.csv
```

The exported CSV has 3,233 rows before the scoring module's annotation filter. Expected SHA-256: `098d8c21cf43d486a548dfc6e4ecdca89221254d6b378ab331973efe98ad4905`. Retain the [upstream license](https://github.com/jinworks/CellChat/blob/592ba49cf4b72130e27eb47a217faf5f84a54f07/LICENSE) and resource citation. The directional table supplies interaction definitions for signaling scores.

## 4. Set the three input paths

Use absolute paths, substituting your own metadata location:

```bash
export SCKNIFE_COUNTS="$PWD/tutorial/inputs/bassez_cohort2_counts.h5ad"
export SCKNIFE_METADATA="/absolute/path/1871-BIOKEY_metaData_cohort2_web.csv"
export SCKNIFE_INTERACTION_DB="$PWD/tutorial/inputs/interactions.csv"
```

Keep the same terminal and activated environment for the {doc}`quickstart`. The scripts check and align cell identifiers; metadata must match the counts. Patient IDs and `Pre`/`On` timepoints are retained for grouping and interpretation. E/NE labels describe clonotype expansion, not clinical response.

The cohort inputs stay on your machine. The metabolic example downloads two small, pinned scCellFie task-support tables to `task_support/` in the output directory and verifies their checksums. For offline use, place `Task_by_Rxn.csv` and `Task_by_Gene.csv` from the [pinned human task resources](https://github.com/earmingol/scCellFie/tree/c66c6be030a8ac871d998f17e592477b97781bf9/task_data/homo_sapiens) there first. See {doc}`graph-resources` for the resource license and citation.
