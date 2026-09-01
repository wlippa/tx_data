"""Per-table build modules: source (TSV/etc.) → canonical Parquet."""

from tx_data.builds import alphamissense, clinical, driver_list, muttable, wgd_calls

__all__ = ["muttable", "wgd_calls", "alphamissense", "clinical", "driver_list"]

BUILDERS = {
    "muttable": muttable.build,
    "wgd_calls": wgd_calls.build,
    "alphamissense": alphamissense.build,
    "clinical": clinical.build,
    "driver_list": driver_list.build,
}
