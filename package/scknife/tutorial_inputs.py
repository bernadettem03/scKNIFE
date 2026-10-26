"""Prepare external inputs for the documented tutorial."""
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.request

CELLCHAT_URL = (
    "https://raw.githubusercontent.com/jinworks/CellChat/"
    "592ba49cf4b72130e27eb47a217faf5f84a54f07/data/CellChatDB.human.rda"
)
INTERACTION_SHA256 = "098d8c21cf43d486a548dfc6e4ecdca89221254d6b378ab331973efe98ad4905"


def prepare_signaling_input(directory, input_path=None):
    """Return a supplied interaction table, or download and cache the tutorial table.

    The automatic export requires base R on its first run. An explicit input_path
    is used as supplied, allowing custom interaction definitions and offline use.
    Existing default files must match the tutorial table and are never overwritten.
    """
    if input_path:
        path = Path(input_path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Signaling input not found: {path}")
        return path

    directory = Path(directory).expanduser().resolve()
    path = directory / "interactions.csv"
    if path.exists():
        if path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == INTERACTION_SHA256:
            return path
        raise ValueError(
            f"{path} differs from the tutorial input. Set SCKNIFE_INTERACTION_DB "
            "to use your own table, or choose an empty input directory."
        )
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError(
            "Install R and make Rscript available to prepare the signaling input, "
            "or set SCKNIFE_INTERACTION_DB to an existing interactions.csv file."
        )
    directory.mkdir(parents=True, exist_ok=True)
    print("Preparing signaling input...")
    with tempfile.TemporaryDirectory(prefix=".signaling-", dir=directory) as work:
        source = Path(work) / "CellChatDB.human.rda"
        output = Path(work) / "interactions.csv"
        with urllib.request.urlopen(CELLCHAT_URL, timeout=120) as response, source.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        code = '''args <- commandArgs(trailingOnly=TRUE)
e <- new.env()
load(args[1], envir=e)
x <- e$CellChatDB.human$interaction
stopifnot(nrow(x)==3233, all(c('ligand.symbol','receptor.symbol','annotation') %in% names(x)))
write.csv(x, args[2])
'''
        subprocess.run(
            [rscript, "--vanilla", "-e", code, str(source), str(output)],
            check=True, capture_output=True, text=True, timeout=120,
        )
        if hashlib.sha256(output.read_bytes()).hexdigest() != INTERACTION_SHA256:
            raise ValueError("The prepared signaling input did not match the expected tutorial table.")
        # Hard-link publication is atomic and does not overwrite a concurrently created file.
        try:
            os.link(output, path)
        except FileExistsError:
            if hashlib.sha256(path.read_bytes()).hexdigest() != INTERACTION_SHA256:
                raise ValueError(f"A different signaling input already exists at {path}")
    return path
