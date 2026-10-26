"""Execute every tutorial cell in the active Python environment."""
import os
import argparse
from pathlib import Path
import tempfile
import nbformat
from nbclient import NotebookClient
from ipykernel.kernelspec import install

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--regression", action="store_true", help="Run the small bundled test example")
args = parser.parse_args()
name = "scKNIFE_regression_1000.ipynb" if args.regression else "scKNIFE_Bassez_Cohort_B_complete_tutorial.ipynb"
source = root / "tutorial" / name
out = Path(os.environ.get("SCKNIFE_TUTORIAL_OUTPUT", root / ("tutorial/outputs" if args.regression else "tutorial/outputs_10000"))).resolve()
out.mkdir(parents=True, exist_ok=True)
os.environ["SCKNIFE_TUTORIAL_OUTPUT"] = str(out)
nb = nbformat.read(source, as_version=4)
# A temporary kernelspec prevents accidentally using a different environment.
with tempfile.TemporaryDirectory(prefix="scknife-kernel-") as tmp:
    install(prefix=tmp, kernel_name="scknife-tutorial")
    os.environ["JUPYTER_PATH"] = str(Path(tmp) / "share/jupyter") + os.pathsep + os.environ.get("JUPYTER_PATH", "")
    client = NotebookClient(nb, timeout=1800, kernel_name="scknife-tutorial",
                            resources={"metadata": {"path": str(root)}}, allow_errors=False)
    client.execute()
nbformat.write(nb, out / name)
print(f"All tutorial cells passed. Outputs: {out}")
