"""Load YAML catalog entries."""

from __future__ import annotations

from pathlib import Path

import yaml

from tx_data.paths import REPO_ROOT

CATALOG_DIR = REPO_ROOT / "catalog"


def load(name: str) -> dict:
    """Load a single catalog entry by table name (no .yml extension)."""
    path = CATALOG_DIR / f"{name}.yml"
    with path.open() as f:
        return yaml.safe_load(f)


def load_all() -> dict[str, dict]:
    """Load every catalog entry keyed by table name."""
    out: dict[str, dict] = {}
    for path in sorted(CATALOG_DIR.glob("*.yml")):
        if path.name.startswith("_"):
            continue
        with path.open() as f:
            out[path.stem] = yaml.safe_load(f)
    return out
