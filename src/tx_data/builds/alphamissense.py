"""Build alphamissense → canonical Parquet.

Reads per-tumour VEP output files (with `##` comment header), skips comments,
strips the leading `#` from the header row, adds `tumour_id` (normalised to
canonical `LTX0001-Tumour1` form), parses `Location` into `chr` + `pos`,
dedupes to per-(tumour, chr, pos, alt) with a single am_pathogenicity per
group (all rows within a group carry the same score by construction —
confirmed 2026-08-31).

File layout (updated 2026-09-02):
  <root>/output/annotated_muttables/tx842/<tumour_id>_muttable_alpha.csv
where `<tumour_id>` in the filename is the muttable spelling
(`LTX0001_tumour1`), which the loader normalises back to canonical form.
"""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl

from tx_data.builds._base import canonical_output_path, log
from tx_data.normalize import canonical_tumour_id
from tx_data.sources import resolve_source, source_entry

TABLE = "alphamissense"


def _read_vep(p: Path) -> pl.DataFrame:
    """Read one VEP file: skip `##` lines, strip `#` from header, load as TSV.

    Files may be named `.tsv` or `.csv`; the body is tab-separated regardless
    (VEP `--tab` output).
    """
    text = p.read_text()
    lines = [ln for ln in text.splitlines() if not ln.startswith("##")]
    if lines and lines[0].startswith("#"):
        lines[0] = lines[0].lstrip("#")
    body = "\n".join(lines)
    return pl.read_csv(
        source=body.encode(),
        separator="\t",
        null_values=["-"],
        infer_schema_length=10000,
    )


def _discover_am_files() -> list[tuple[str, Path]]:
    """Return ``[(tumour_id_from_filename, path), …]`` for every AM file on disk.

    Uses ``relative_path_pattern`` from ``_sources.yml`` to derive both the
    directory to scan and the filename shape to parse — so a future rename
    only needs the yml change.
    """
    entry = source_entry(TABLE)
    pattern = entry["relative_path_pattern"]  # e.g. output/annotated_muttables/tx842/{tumour_id}_muttable_alpha.csv
    probe = resolve_source(TABLE, tumour_id="__probe__")
    scan_dir = probe.parent
    filename_template = Path(pattern).name  # e.g. `{tumour_id}_muttable_alpha.csv`
    prefix, suffix = filename_template.split("{tumour_id}", 1)
    fname_re = re.compile(re.escape(prefix) + r"(?P<tid>.+?)" + re.escape(suffix) + r"$")

    if not scan_dir.is_dir():
        log(f"AM scan dir does not exist: {scan_dir}")
        return []

    found: list[tuple[str, Path]] = []
    for p in sorted(scan_dir.iterdir()):
        if not p.is_file():
            continue
        m = fname_re.match(p.name)
        if not m:
            continue
        found.append((m.group("tid"), p))
    return found


def build() -> pl.DataFrame:
    discovered = _discover_am_files()
    log(f"discovered {len(discovered)} AM files")

    frames = []
    for tid_from_filename, p in discovered:
        tid_canonical = canonical_tumour_id(tid_from_filename)
        df = _read_vep(p).with_columns(pl.lit(tid_canonical).alias("tumour_id"))
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
