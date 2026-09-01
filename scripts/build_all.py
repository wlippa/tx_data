"""Build every canonical Parquet from source, in dependency order."""
from tx_data.builds import BUILDERS

if __name__ == "__main__":
    for name in ("muttable", "clinical", "wgd_calls", "alphamissense", "driver_list"):
        BUILDERS[name]()
