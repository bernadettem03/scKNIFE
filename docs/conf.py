"""Sphinx configuration for the scKNIFE user documentation."""

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "package"))

project = "scKNIFE"
author = "Bernadette Mathew"
copyright = "2026, Bernadette Mathew"
try:
    release = version("scknife")
except PackageNotFoundError:
    release = "0.1.0"
version = release

extensions = [
    "myst_parser", "nbsphinx", "nbsphinx_link", "sphinx.ext.autodoc",
    "sphinx.ext.autosummary", "sphinx.ext.napoleon",
    "sphinx.ext.viewcode", "sphinx_copybutton", "sphinx_design",
]
autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = True
myst_enable_extensions = ["colon_fence", "deflist", "fieldlist", "substitution"]
nbsphinx_execute = "never"
nbsphinx_allow_errors = False
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "downloads/*.ipynb", "licenses/*.md", "tutorials/bassez-cohort-b.nblink"]
suppress_warnings = ["config.cache"]

html_theme = "furo"
html_title = "scKNIFE documentation"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
# Set an initial site preference before Furo reads the browser's theme setting.
html_js_files = [("light-default.js", {"priority": 100})]
pygments_style = "friendly"
html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "source_repository": "https://github.com/bernadettem03/scKNIFE/",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#126b74",
        "color-brand-content": "#126b74",
        "color-background-primary": "#ffffff",
        "color-background-secondary": "#f7f9fb",
        "color-sidebar-background": "#f7f9fb",
        "color-foreground-primary": "#243344",
        "color-foreground-secondary": "#526171",
    },
}
