"""Render executable tutorial pages into a notebook (preserve matching outputs)."""

from pathlib import Path
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "tutorial/scKNIFE_Bassez_Cohort_B_complete_tutorial.ipynb"
PAGES = [
    "setup",
    "clustering",
    "factors",
    "metabolic",
    "go",
    "signaling",
    "states",
    "validation",
]


def make_notebook(pages=PAGES):
    cells = []

    def add(kind, source, stage):
        if not source.strip():
            return
        cell = dict(
            cell_type=kind,
            source=source.strip().splitlines(keepends=True),
            metadata={"tutorial_stage": stage},
        )
        if kind == "code":
            cell.update(outputs=[], execution_count=None)
        cells.append(cell)

    add("markdown", "# scKNIFE: 36,721-cell tutorial", "installation")
    add(
        "markdown",
        '## Installation\n\nUse Python 3.9.21. From the cloned repository, install the tutorial environment:\n\n```bash\npython -m pip install -e ".[tutorial]"\n```',
        "installation",
    )
    add(
        "code",
        '%matplotlib inline\nimport scknife\nfrom importlib.metadata import version\nprint(version("scknife"))',
        "installation",
    )
    for stage in pages:
        page = (ROOT / f"docs/tutorials/36721-{stage}.md").read_text()
        page = re.sub(r"^\{download\}.*$", "", page, flags=re.M)
        page = re.sub(r"^\([\w-]+\)=\s*$", "", page, flags=re.M)
        page = re.sub(r"^---\n.*?---\n", "", page, flags=re.S)
        page = re.sub(r"```\{(?:image|csv-table)\}.*?```", "", page, flags=re.S)
        page = re.sub(
            r"\]\(([^(): ]+\.md)(#[^)]*)?\)",
            lambda m: "](" + ("../docs/tutorials/" + m[1]) + (m[2] or "") + ")",
            page,
        )
        if stage == "factors":
            add(
                "markdown",
                "## Analysis tasks\n\nChoose any task using the saved model and activities.",
                stage,
            )
        if stage == "states":
            add("markdown", "## Optional examples", stage)
        for chunk in re.split(r"(```python\n.*?```)", page, flags=re.S):
            if chunk.startswith("```python\n"):
                add("code", chunk[10:-3], stage)
            else:
                add(
                    "markdown",
                    re.sub(
                        r"^(#+) ",
                        lambda m: (
                            ("#" if stage in ("setup", "clustering") else "##")
                            + m[1]
                            + " "
                        ),
                        chunk,
                        flags=re.M,
                    ),
                    stage,
                )
    return dict(
        cells=cells,
        metadata={
            "kernelspec": {
                "name": "python3",
                "display_name": "Python 3 (scKNIFE)",
                "language": "python",
            },
            "language_info": {"name": "python", "version": "3.9.21"},
        },
        nbformat=4,
        nbformat_minor=4,
    )


if __name__ == "__main__":
    nb = make_notebook()
    if "--check" in sys.argv:
        old = json.loads(TARGET.read_text())
        assert [(c["cell_type"], c["source"]) for c in old["cells"]] == [
            (c["cell_type"], c["source"]) for c in nb["cells"]
        ]
    else:
        TARGET.write_text(json.dumps(nb, indent=1) + "\n")
