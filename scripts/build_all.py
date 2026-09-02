"""Build every canonical Parquet from source, in dependency order.

Writes a `data/_build_manifest.json` at the end recording the data_root
used, the TX_DATA_ROOT env var seen, whether this was a mock build, the
timestamp, and per-table row counts. Downstream consumers (e.g.
mut_essential_wgd/run_analysis.py) read this to prove which cohort a
given set of parquets came from.
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import polars as pl

from tx_data.builds import BUILDERS
from tx_data.builds._base import DATA_DIR, log
from tx_data.paths import data_root, data_root_is_mock

TABLES = ("muttable", "clinical", "wgd_calls", "alphamissense", "driver_list")


def _write_manifest() -> None:
    row_counts: dict[str, int] = {}
    for name in TABLES:
        p = DATA_DIR / f"{name}.parquet"
        if p.is_file():
            row_counts[name] = pl.read_parquet(p).height
        # clinical builder additionally emits clinical_per_tumour; capture it.
        if name == "clinical":
            p2 = DATA_DIR / "clinical_per_tumour.parquet"
            if p2.is_file():
                row_counts["clinical_per_tumour"] = pl.read_parquet(p2).height

    manifest = {
        "built_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "data_root": str(data_root()),
        "data_root_is_mock": data_root_is_mock(),
        "env_TX_DATA_ROOT": os.environ.get("TX_DATA_ROOT"),
        "table_row_counts": row_counts,
    }
    manifest_path = DATA_DIR / "_build_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    log(f"wrote {manifest_path}")
    log(
        f"data_root={manifest['data_root']} "
        f"is_mock={manifest['data_root_is_mock']} "
        f"tables={row_counts}"
    )


if __name__ == "__main__":
    for name in TABLES:
        BUILDERS[name]()
    _write_manifest()
