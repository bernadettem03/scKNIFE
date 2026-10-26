# Contributing

Use GitHub issues to report problems or suggest improvements.

## Build the documentation

Install scKNIFE in the Python environment described in the [installation guide](docs/installation.md), then run:

```bash
pip install -r docs/requirements.txt
python scripts/check_tutorial_contract.py
sphinx-build -W --keep-going -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` to preview the documentation.
