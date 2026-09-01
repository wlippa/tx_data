"""End-to-end smoke test: mock → build → assert invariants on canonical Parquets."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "src"
SCRIPTS = REPO / "scripts"
DATA = REPO / "data"


def _run(script: Path) -> None:
    r = subprocess.run(
        [sys.executable, str(script)],
        cwd=REPO,
        env={"PYTHONPATH": str(SRC), "PATH": "/usr/bin:/bin"},
        check=True,
        capture_output=True,
    )
    assert r.returncode == 0, r.stderr.decode()


@pytest.fixture(scope="module", autouse=True)
def _build_all():
    _run(SCRIPTS / "generate_mocks.py")
    _run(SCRIPTS / "build_all.py")


def test_muttable_shape():
    df = pl.read_parquet(DATA / "muttable.parquet")
    # Canonical tumour_id form applied.
    assert df["patient_tumour"].str.contains("-Tumour").all()
    # Clone form: string 'cloneN' or null.
    assert df["mutation_cluster"].drop_nulls().str.starts_with("clone").all()


def test_wgd_calls_derived_columns():
    df = pl.read_parquet(DATA / "wgd_calls.parquet")
    for col in ("is_root", "is_wgd_event_clone",
                "cumulative_gds_to_clone", "has_wgd_ancestor"):
        assert col in df.columns
    # Cumulative WGD at MRCA equals number_of_gds_at_clone at MRCA.
    roots = df.filter(pl.col("is_root"))
    assert (roots["cumulative_gds_to_clone"] == roots["number_of_gds_at_clone"]).all()
    # Cumulative WGD is non-decreasing along parent→child (spot-check on one tumour).
    for tumour_id, sub in df.group_by("tumour_id"):
        cmap = dict(zip(sub["clone"], sub["cumulative_gds_to_clone"]))
        pmap = dict(zip(sub["clone"], sub["parent"]))
        for c, p in pmap.items():
            if p == "diploid":
                continue
            assert cmap[c] >= cmap[p], f"{tumour_id[0]}: {c}({cmap[c]}) < {p}({cmap[p]})"


def test_wgd_calls_class_values():
    df = pl.read_parquet(DATA / "wgd_calls.parquet")
    resolved = df.filter(pl.col("status") == "resolved")
    seen_classes = set(resolved["class"].drop_nulls().unique().to_list())
    # Mock cohort covers all four class values.
    assert {"no_wgd", "clonal_wgd", "mut_supported", "ploidy_only"}.issubset(seen_classes)


def test_alphamissense_dedup():
    df = pl.read_parquet(DATA / "alphamissense.parquet")
    # One row per (tumour_id, chr, pos, alt).
    n = df.select(["tumour_id", "chr", "pos", "alt"]).n_unique()
    assert n == df.height


def test_clinical_per_tumour_unpivot():
    per_p = pl.read_parquet(DATA / "clinical.parquet")
    per_t = pl.read_parquet(DATA / "clinical_per_tumour.parquet")
    # Multi-tumour patient LTX0008 should have 2 rows in per_tumour.
    assert per_p.filter(pl.col("Patient_ID") == "LTX0008").height == 1
    assert per_t.filter(pl.col("Patient_ID") == "LTX0008").height == 2
    # tumour_id in canonical form.
    assert per_t["tumour_id"].str.contains("-Tumour").all()


def test_id_normalisation_consistency():
    """The tumour_id form in muttable must match wgd_calls after normalisation."""
    mut = pl.read_parquet(DATA / "muttable.parquet")
    wgd = pl.read_parquet(DATA / "wgd_calls.parquet")
    mut_tumours = set(mut["patient_tumour"].unique().to_list())
    wgd_tumours = set(wgd["tumour_id"].unique().to_list())
    # muttable side already normalised by builder.
    common = mut_tumours & wgd_tumours
    assert common, f"no overlap; mut={list(mut_tumours)[:3]} wgd={list(wgd_tumours)[:3]}"
