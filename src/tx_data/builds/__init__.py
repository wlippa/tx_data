"""Per-table build modules: source (TSV/etc.) → canonical Parquet."""

from tx_data.builds import alphamissense, clinical, muttable, wgd_calls

__all__ = ["muttable", "wgd_calls", "alphamissense", "clinical"]

BUILDERS = {
    "muttable": muttable.build,
    "wgd_calls": wgd_calls.build,
    "alphamissense": alphamissense.build,
    "clinical": clinical.build,
}
