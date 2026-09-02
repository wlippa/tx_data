"""Build muttable → canonical Parquet.

Applies:
  - tumour_id normalisation (`LTX0001_tumour1` → `LTX0001-Tumour1`)
  - clone normalisation on `mutation_cluster` (float → `clone{N}` string)
  - dtype coercion per the catalog

Source: TSV gz. Output: single Parquet at `data/muttable.parquet`.
"""

from __future__ import annotations

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.normalize import canonical_clone_expr, canonical_tumour_id_expr
from tx_data.sources import resolve_source

TABLE = "muttable"


def build() -> pl.DataFrame:
    src = resolve_source(TABLE)
    log(f"read {src}")

    # Force `chr` to String at READ time. If we let polars infer, the first
    # 10 000 rows may be all autosomes → i64 inferred → parse fails when it
    # hits X / Y / MT further into the file.
    df = pl.read_csv(
        src,
        separator="\t",
        null_values=["NA", ""],
        infer_schema_length=10000,
        schema_overrides={"chr": pl.String},
    )

    df = df.with_columns(
        [
            # Canonicalise tumour_id from muttable's `patient_tumour` form.
            canonical_tumour_id_expr("patient_tumour"),
            # Canonicalise clone (was float mutation_cluster).
            canonical_clone_expr("mutation_cluster"),
        ]
    )

    out = canonical_output_path(TABLE)
    df.write_parquet(out)
    log(f"wrote {out} ({df.height} rows × {df.width} cols)")
    return df


if __name__ == "__main__":
    build()
