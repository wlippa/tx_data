"""Resolve the data root. TX_DATA_ROOT env var wins; else config/paths.yml."""

from __future__ import annotations

import os
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = REPO_ROOT / "config" / "paths.yml"


def data_root() -> Path:
    env = os.environ.get("TX_DATA_ROOT")
    if env:
        return Path(env)
    with CONFIG_PATH.open() as f:
        cfg = yaml.safe_load(f)
    root = Path(cfg["data_root"])
    if not root.is_absolute():
        root = REPO_ROOT / root
    return root


def derived_dir() -> Path:
    return REPO_ROOT / "data"


def db_path() -> Path:
    return REPO_ROOT / "db" / "tx_data.duckdb"
