"""Capture long training logs without hiding exceptions."""

from contextlib import redirect_stdout
from functools import wraps


def _run_log(relative_path):
    """Keep iterative fit output in a local file; errors still propagate."""

    def decorate(function):
        @wraps(function)
        def wrapped(*args, **kwargs):
            from .core import output_dir

            path = output_dir() / relative_path
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("w") as log, redirect_stdout(log):
                return function(*args, **kwargs)

        return wrapped

    return decorate
