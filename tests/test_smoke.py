"""Smoke tests: path resolution + catalog loading."""

import os
from pathlib import Path

from tx_data import catalog, paths


def test_repo_root_exists():
    assert paths.REPO_ROOT.exists()


def test_data_root_returns_path():
    root = paths.data_root()
    assert isinstance(root, Path)


def test_data_root_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("TX_DATA_ROOT", str(tmp_path))
    assert paths.data_root() == tmp_path


def test_data_root_default_is_nemo_mock(monkeypatch):
    monkeypatch.delenv("TX_DATA_ROOT", raising=False)
    root = paths.data_root()
    assert root.name == "nemo_mock"


def test_catalog_dir_exists():
    assert catalog.CATALOG_DIR.exists()


def test_catalog_load_all_returns_dict():
    assert isinstance(catalog.load_all(), dict)
