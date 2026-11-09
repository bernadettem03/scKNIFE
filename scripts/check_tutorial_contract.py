"""Check tutorial sources, executable examples, notebook and graph references."""

import ast
import json
import re
from pathlib import Path
from build_main_notebook import ROOT, TARGET, PAGES, make_notebook

suffixes = [
    "",
    "_organization",
    "_factors",
    "_metabolic",
    "_go",
    "_signaling",
    "_states",
    "_validation",
]
download = {
    "metabolic": "metabolic_walkthrough",
    "go": "go_differential",
    "factors": "factor_interpretation",
    "signaling": "signaling_walkthrough",
}
for stage, suffix in zip(PAGES, suffixes):
    page = ROOT / f"docs/tutorials/36721-{stage}.md"
    text = page.read_text()
    code = (
        "\n\n".join(i.rstrip() for i in re.findall(r"```python\n(.*?)```", text, re.S))
        + "\n"
    )
    assert not any(
        isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
        for n in ast.walk(ast.parse(code))
    ), stage
    assert code == (ROOT / f"scripts/run_tutorial_36721{suffix}.py").read_text(), stage
    for image in re.findall(r"```\{image\} (.+)", text):
        assert (page.parent / image).is_file(), image
    if stage in download:
        assert code == (ROOT / f"docs/downloads/{download[stage]}.py").read_text()
expected = make_notebook()
nb = json.loads(TARGET.read_text())
assert [(c["cell_type"], c["source"]) for c in nb["cells"]] == [
    (c["cell_type"], c["source"]) for c in expected["cells"]
]
assert (
    TARGET.read_bytes()
    == (ROOT / "docs/downloads/scKNIFE_36721_tutorial.ipynb").read_bytes()
)
print("All tutorial stages, downloads and notebook sources agree.")
