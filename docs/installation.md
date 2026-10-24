# Installation and inputs

scKNIFE currently supports **Python 3.9.21** with a reproducibility-focused dependency set. Use a clean environment rather than installing into an existing Scanpy environment.

## Install from GitHub

```bash
conda create -n scknife python=3.9.21 pip
conda activate scknife
pip install "git+https://github.com/bernadettem03/scKNIFE.git"
```

The repository is public; no GitHub authentication is required.

## Editable installation

```bash
git clone https://github.com/bernadettem03/scKNIFE.git
cd scKNIFE
conda create -n scknife python=3.9.21 pip
conda activate scknife
pip install -e .
```

## Run the 10,000-cell tutorial

After cloning the repository, install the tutorial dependencies:

```bash
pip install -e '.[tutorial]'
```

Follow {doc}`tutorial-inputs` to obtain and prepare the external data. Then run the commands in {doc}`quickstart`. The {doc}`tutorial-10000` explains the same analysis in detail. The supplied graphs are loaded unchanged; this workflow does not require graph reconstruction or the graph-building extras.

## Verify

```python
import scknife
from scknife import run_full_pipeline, run_tuning
print(scknife.__file__)
print(run_full_pipeline)
```

## Clustering dependencies

H-space clustering uses Scanpy's Leiden implementation. The `tutorial` extra above installs these dependencies; install them separately only if you installed the base package:

```bash
pip install python-igraph leidenalg
```

## Expected resources

For the documented 10,000-cell example, sampling, tuning and fitting took about 8 minutes and used approximately 3.6 GiB peak process memory on the tested machine. Input conversion and installation are additional. Runtime depends on hardware, graph size and settings.

## Build the documentation locally

```bash
pip install -r docs/requirements.txt
sphinx-build -W --keep-going -b html docs docs/_build/html
```

See {doc}`troubleshooting` for common setup problems.

## Continue the tutorial

Obtain the files in {doc}`tutorial-inputs`, then continue to {doc}`tutorials/10000-setup`. See {doc}`tutorial-10000` for the full numbered sequence.
