# Input files

The tutorial uses Bassez cohort-2 counts and metadata, plus the supplied Version 1 graph. It selects 10,000 TNBC cells automatically.

## Download counts and metadata

Download both files from the **Data files** section of the [Lambrechts lab Bassez data page](https://lambrechtslab.sites.vib.be/en/single-cell):

| File | Contents |
|---|---|
| `1867-counts_cells_cohort2.rds` | Cohort-2 counts |
| `1871-BIOKEY_metaData_cohort2_web.csv` | Matching cell metadata |

## Convert the counts

Use R with the Matrix package and the Python environment from {doc}`installation`:

```bash
mkdir -p tutorial/inputs/converted_counts
Rscript scripts/prepare_bassez_counts.R /absolute/path/1867-counts_cells_cohort2.rds tutorial/inputs/converted_counts
python scripts/prepare_bassez_counts.py tutorial/inputs/converted_counts tutorial/inputs/bassez_cohort2_counts.h5ad
```

Skip this step if you already have the equivalent unnormalized H5AD file.

## Set input paths

Replace the metadata path with your file location:

```bash
export SCKNIFE_COUNTS="$PWD/tutorial/inputs/bassez_cohort2_counts.h5ad"
export SCKNIFE_METADATA="/absolute/path/1871-BIOKEY_metaData_cohort2_web.csv"
```

Continue with {doc}`quickstart` or {doc}`tutorials/10000-setup`.
