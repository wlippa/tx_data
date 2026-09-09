"""Per-table build modules: source (TSV/etc.) → canonical Parquet."""

from tx_data.builds import (
    alpaca,
    alphamissense,
    clinical,
    clone_proportions,
    driver_list,
    kallisto,
    muttable,
    wgd_calls,
)

__all__ = [
    "muttable",
    "wgd_calls",
    "alphamissense",
    "clinical",
    "driver_list",
    "alpaca",
    "kallisto",
    "clone_proportions",
]

BUILDERS = {
    "muttable": muttable.build,
    "wgd_calls": wgd_calls.build,
    "alphamissense": alphamissense.build,
    "clinical": clinical.build,
    "driver_list": driver_list.build,
    "alpaca": alpaca.build,
    "kallisto": kallisto.build,
    "clone_proportions": clone_proportions.build,
}
