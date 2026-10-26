"""Check automatic tutorial inputs without network access or an R installation."""
import hashlib
import io
from pathlib import Path
import subprocess

import pytest
from scknife import tutorial_inputs as inputs


def no_download(*args, **kwargs):
    raise AssertionError("An existing or explicit input must not be downloaded")


def test_explicit_input_is_used_without_download(tmp_path, monkeypatch):
    source = tmp_path / "custom.csv"
    source.write_text("ligand.symbol,receptor.symbol\nA,B\n")
    monkeypatch.setattr(inputs.urllib.request, "urlopen", no_download)
    assert inputs.prepare_signaling_input(tmp_path / "cache", source) == source
    with pytest.raises(FileNotFoundError):
        inputs.prepare_signaling_input(tmp_path / "cache", tmp_path / "missing.csv")


def test_cached_input_works_offline_and_without_r(tmp_path, monkeypatch):
    source = tmp_path / "interactions.csv"
    source.write_bytes(b"verified input")
    monkeypatch.setattr(inputs, "INTERACTION_SHA256", hashlib.sha256(source.read_bytes()).hexdigest())
    monkeypatch.setattr(inputs.urllib.request, "urlopen", no_download)
    monkeypatch.setattr(inputs.shutil, "which", lambda _: None)
    assert inputs.prepare_signaling_input(tmp_path) == source


def test_existing_mismatched_file_is_not_overwritten(tmp_path, monkeypatch):
    source = tmp_path / "interactions.csv"
    source.write_bytes(b"user file")
    monkeypatch.setattr(inputs.urllib.request, "urlopen", no_download)
    with pytest.raises(ValueError, match="differs"):
        inputs.prepare_signaling_input(tmp_path)
    assert source.read_bytes() == b"user file"


def test_first_run_requires_r(tmp_path, monkeypatch):
    monkeypatch.setattr(inputs.shutil, "which", lambda _: None)
    monkeypatch.setattr(inputs.urllib.request, "urlopen", no_download)
    with pytest.raises(RuntimeError, match="Install R"):
        inputs.prepare_signaling_input(tmp_path)


@pytest.mark.parametrize("outcome", ["valid", "invalid", "failed"])
def test_only_verified_complete_exports_are_cached(tmp_path, monkeypatch, outcome):
    content = b"expected interaction export\n"
    monkeypatch.setattr(inputs, "INTERACTION_SHA256", hashlib.sha256(content).hexdigest())
    monkeypatch.setattr(inputs.shutil, "which", lambda _: "/test/Rscript")
    monkeypatch.setattr(inputs.urllib.request, "urlopen", lambda *a, **k: io.BytesIO(b"rdata"))
    def export(command, **kwargs):
        Path(command[-1]).write_bytes(content if outcome == "valid" else b"partial")
        if outcome == "failed":
            raise subprocess.CalledProcessError(1, command)
    monkeypatch.setattr(inputs.subprocess, "run", export)
    if outcome == "valid":
        path = inputs.prepare_signaling_input(tmp_path)
        assert path.read_bytes() == content
    else:
        with pytest.raises((ValueError, subprocess.CalledProcessError)):
            inputs.prepare_signaling_input(tmp_path)
        assert not (tmp_path / "interactions.csv").exists()
    assert not list(tmp_path.glob(".signaling-*"))
