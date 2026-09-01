"""Build alphamissense → canonical Parquet.

Reads per-tumour VEP output files (with `##` comment header), skips comments,
strips the leading `#` from the header row, adds `tumour_id`, parses
`Location` into `chr` + `pos`, dedupes to per-(tumour, chr, pos, alt) with a
single am_pathogenicity per group (all rows within a group carry the same
score by construction — confirmed 2026-08-31).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.paths import data_root
from tx_data.sources import resolve_source, source_entry

TABLE = "alphamissense"


def _read_vep(p: Path) -> pl.DataFrame:
    """Read one VEP file: skip `##` lines, strip `#` from header, load TSV."""
    text = p.read_text()
    lines = [ln for ln in text.splitlines() if not ln.startswith("##")]
    # First surviving line is the header — strip leading `#`.
    if lines and lines[0].startswith("#"):
        lines[0] = lines[0].lstrip("#")
    body = "\n".join(lines)
    return pl.read_csv(
        source=body.encode(),
        separator="\t",
        null_values=["-"],
        infer_schema_length=10000,
    )


def _list_tumour_dirs() -> list[str]:
    """Discover per-tumour subdirs under the alphamissense root."""
    entry = source_entry(TABLE)
    # Enumerate tumour_ids by listing the parent of the pattern.
    # relative_path_pattern is `annotated_muttables/tx842/{tumour_id}/{tumour_id}_muttable_annotated.tsv`.
    pattern = entry["relative_path_pattern"]
    parent_rel = pattern.split("{tumour_id}")[0]  # `annotated_muttables/tx842/`
    src_root = resolve_source(TABLE, tumour_id="__probe__").parent.parent  # up to `tx842/`
    if not src_root.is_dir():
        # In mock, the pattern's parent may be at data_root/alt/... — fall back.
        src_root = (data_root() / "alt" / entry.get("source_root", "") / parent_rel).resolve()
    return sorted(p.name for p in src_root.iterdir() if p.is_dir())


def build() -> pl.DataFrame:
    tumours = _list_tumour_dirs()
    log(f"discovered {len(tumours)} tumour dirs")

    frames = []
    for t in tumours:
        p = resolve_source(TABLE, tumour_id=t)
        if not p.is_file():
            log(f"  skip {t}: no AM file at {p}")
            continue
        df = _read_vep(p).with_columns(pl.lit(t).alias("tumour_id"))
        frames.append(df)

    raw = pl.concat(frames, how="vertical_relaxed") if frames else pl.DataFrame()

    if raw.is_empty():
        log("no AM rows to write")
        out = canonical_output_path(TABLE)
        raw.write_parquet(out)
        return raw

    # Parse Location "chr:pos" or "chr:start-stop".
    loc_split = pl.col("Location").str.split_exact(":", 1)
    parsed = raw.with_columns(
        [
            loc_split.struct.field("field_0").alias("chr"),
            loc_split.struct.field("field_1")
                .str.split_exact("-", 1)
                .struct.field("field_0")
                .cast(pl.Int64)
                .alias("pos"),
            pl.col("Allele").alias("alt"),
            pl.col("am_pathogenicity").cast(pl.Float64, strict=False),
        ]
    )

    # Dedup to per-variant AM score. `max` collapses agreeing values and ignores nulls.
    scored = parsed.group_by(["tumour_id", "chr", "pos", "alt"]).agg(
        [
            pl.col("am_pathogenicity").max().alias("am_pathogenicity"),
            pl.col("am_class").drop_nulls().first().alias("am_class"),
            pl.col("Consequence").drop_nulls().first().alias("consequence_any"),
            pl.col("Gene").drop_nulls().first().alias("gene_ensembl_any"),
        ]
    )

    out = canonical_output_path(TABLE)
    scored.write_parquet(out)
    log(f"wrote {out} ({scored.height} rows × {scored.width} cols)")
    return scored


if __name__ == "__main__":
    build()
