"""Generate synthetic mock cohort into ``nemo_mock/``.

Produces four coherent tables spanning every meaningful edge case:

  - 10 tumours (~9 patients, one with a second tumour for multi-tumour coverage)
  - Every wgd_calls `class`: no_wgd, clonal_wgd, mut_supported, ploidy_only
  - A `sub_class` split (parallel / sequential / single subclonal GD)
  - Every status: mostly `resolved`, one `unresolved`, one `needs_follow_up`
  - Both High and Low quality
  - A subclonal-WGD tumour where the WGD event is on an internal clone (so
    tree traversal has to run) plus a parallel-WGD tumour with two independent
    doublings
  - Every gene set: essential-only, driver-only, non-essential (drawn from
    fixed lists so we can assert on them in tests)
  - Every variant class: SNV missense, SNV synonymous, SNV nonsense, indel

Files written under ``nemo_mock/``:
  - ``tx842_mutation_table.tsv.gz``  (long format, one row per mutation×sample)
  - ``WGD/release_tx842/WGD_calls.tsv``
  - ``tx842_clinical.tsv``
  - ``alt/cn_ccf_alphamissense/annotated_muttables/tx842/<tumour>/<tumour>_muttable_annotated.tsv``

Deterministic (fixed seed) so mocks re-generate identically.
"""

from __future__ import annotations

import gzip
import random
from dataclasses import dataclass, field
from pathlib import Path

from tx_data.paths import REPO_ROOT

MOCK_ROOT = REPO_ROOT / "nemo_mock"
AM_ROOT = MOCK_ROOT / "alt" / "cn_ccf_alphamissense" / "annotated_muttables" / "tx842"
WGD_DIR = MOCK_ROOT / "WGD" / "release_tx842"

RNG = random.Random(20260831)

# --- Fixed gene panel (chr:pos ranges chosen to look hg19-ish) ---------------
# Each gene gets a fake 1kb coding region on chr10 / chr17 etc. Keep it tiny so
# mock TSVs stay small.
GENES = [
    # (symbol, chr, start_pos, ensembl_id, category)
    ("TP53",     "17",  7577000, "ENSG00000141510", "driver"),
    ("KRAS",     "12", 25398000, "ENSG00000133703", "driver"),
    ("PDCD11",   "10", 105200000, "ENSG00000148843", "nonessential"),
    ("GSTO2",    "10", 106039000, "ENSG00000065621", "nonessential"),
    ("ARHGAP19", "10",  99003000, "ENSG00000213390", "nonessential"),
    ("EGFR",     "7",  55086000, "ENSG00000146648", "driver"),
    ("MYC",      "8", 128748000, "ENSG00000136997", "driver"),
    ("POLR2A",   "17",  7389000, "ENSG00000181222", "essential"),
    ("RPL5",     "1",  93297000, "ENSG00000122406", "essential"),
    ("EIF3A",    "10", 120790000, "ENSG00000107581", "essential"),
    ("RPS3",     "11", 75111000, "ENSG00000149273", "essential"),
    ("MRPS2",    "9", 139036000, "ENSG00000122140", "essential"),
    ("HSPE1",    "2", 198365000, "ENSG00000115541", "essential"),
    ("SNRPD1",   "18", 19195000, "ENSG00000167088", "essential"),
    ("PSMA1",    "11", 14528000, "ENSG00000129084", "essential"),
]

CODON_TABLE = {
    "missense":   ("Ctt", "Gtt", "L", "V"),
    "synonymous": ("Ctt", "Ctc", "L", "L"),
    "nonsense":   ("Cga", "Tga", "R", "*"),
}


# --- Tumour cohort definition -------------------------------------------------

@dataclass
class Clone:
    name: str
    parent: str
    gds_at_clone: int = 0


@dataclass
class TumourSpec:
    patient_id: str
    tumour_ordinal: int
    histology_group: str
    class_: str
    sub_class: str | None
    status: str
    quality: str
    n_samples: int
    clones: list[Clone]
    total_gds: int = field(init=False)

    def __post_init__(self) -> None:
        self.total_gds = sum(c.gds_at_clone for c in self.clones)

    @property
    def tumour_id_muttable(self) -> str:
        return f"{self.patient_id}_tumour{self.tumour_ordinal}"

    @property
    def tumour_id_canonical(self) -> str:
        return f"{self.patient_id}-Tumour{self.tumour_ordinal}"


def _linear_tree(n_clones: int) -> list[Clone]:
    """clone1 (root) → clone2 → clone3 → ..."""
    clones = [Clone("clone1", "diploid")]
    for i in range(2, n_clones + 1):
        clones.append(Clone(f"clone{i}", f"clone{i - 1}"))
    return clones


def _branching_tree() -> list[Clone]:
    """
        clone1
       /      \\
    clone2    clone3
       |         |
    clone4    clone5
    """
    return [
        Clone("clone1", "diploid"),
        Clone("clone2", "clone1"),
        Clone("clone3", "clone1"),
        Clone("clone4", "clone2"),
        Clone("clone5", "clone3"),
    ]


def cohort() -> list[TumourSpec]:
    specs: list[TumourSpec] = []

    # 1. No-WGD LUAD.
    specs.append(TumourSpec("LTX0001", 1, "LUAD", "no_wgd", None,
                            "resolved", "High", 3, _linear_tree(3)))

    # 2. Clonal-WGD LUAD (single WGD at MRCA).
    tree = _linear_tree(4)
    tree[0].gds_at_clone = 1
    specs.append(TumourSpec("LTX0002", 1, "LUAD", "clonal_wgd", None,
                            "resolved", "High", 3, tree))

    # 3. Clonal-WGD LUSC.
    tree = _linear_tree(3)
    tree[0].gds_at_clone = 1
    specs.append(TumourSpec("LTX0003", 1, "LUSC", "clonal_wgd", None,
                            "resolved", "Low", 2, tree))

    # 4. Subclonal WGD (mut_supported), single event on an internal clone.
    tree = _branching_tree()
    tree[3].gds_at_clone = 1  # clone4 acquires the WGD
    specs.append(TumourSpec("LTX0004", 1, "LUAD", "mut_supported",
                            "single subclonal GD", "resolved", "High", 4, tree))

    # 5. Parallel subclonal WGDs (mut_supported, parallel only).
    tree = _branching_tree()
    tree[3].gds_at_clone = 1  # clone4
    tree[4].gds_at_clone = 1  # clone5 — independent branch
    specs.append(TumourSpec("LTX0005", 1, "LUAD", "mut_supported",
                            "parallel only", "resolved", "High", 4, tree))

    # 6. Sequential subclonal WGDs (mut_supported, sequential only).
    tree = _linear_tree(4)
    tree[1].gds_at_clone = 1  # clone2
    tree[2].gds_at_clone = 1  # clone3 — sequential on the same path
    specs.append(TumourSpec("LTX0006", 1, "LUAD", "mut_supported",
                            "sequential only", "resolved", "High", 3, tree))

    # 7. ploidy_only subclonal WGD (weaker call).
    tree = _branching_tree()
    tree[3].gds_at_clone = 1
    specs.append(TumourSpec("LTX0007", 1, "LUAD", "ploidy_only",
                            "single subclonal GD", "resolved", "Low", 3, tree))

    # 8. Multi-tumour patient: Tumour1 is clonal_wgd, Tumour2 is no_wgd.
    tree = _linear_tree(3)
    tree[0].gds_at_clone = 1
    specs.append(TumourSpec("LTX0008", 1, "LUAD", "clonal_wgd", None,
                            "resolved", "High", 2, tree))
    specs.append(TumourSpec("LTX0008", 2, "LUSC", "no_wgd", None,
                            "resolved", "High", 2, _linear_tree(2)))

    # 9. Unresolved tumour — should be filtered out downstream.
    specs.append(TumourSpec("LTX0009", 1, "other", "no_wgd", None,
                            "unresolved", "Low", 2, _linear_tree(2)))

    # 10. needs_follow_up tumour — also filtered out.
    tree = _linear_tree(3)
    tree[0].gds_at_clone = 2  # double-clonal at MRCA (exercises number_of_gds > 1)
    specs.append(TumourSpec("LTX0010", 1, "LUAD", "clonal_wgd", None,
                            "needs_follow_up", "High", 2, tree))

    return specs


# --- Muttable ----------------------------------------------------------------

def _assign_mutations(spec: TumourSpec) -> list[dict]:
    """Return a list of mutation-clone assignments for this tumour.

    Every non-root clone gets 2–4 mutations across the fixed gene panel.
    Root gets more (trunk mutations).
    """
    muts: list[dict] = []
    for i, clone in enumerate(spec.clones):
        n = 6 if clone.parent == "diploid" else RNG.randint(2, 4)
        for _ in range(n):
            g = RNG.choice(GENES)
            # Pick a variant class.
            klass = RNG.choices(
                ["missense", "synonymous", "nonsense", "indel"],
                weights=[0.55, 0.25, 0.10, 0.10],
            )[0]
            pos = g[2] + RNG.randint(0, 990)
            if klass == "indel":
                # Simple 1bp deletion.
                start = pos
                stop = pos
                pos_left = pos - 1
                ref = RNG.choice(["A", "C", "G", "T"])
                var = "-"
                is_snv, is_indel, is_mnv = False, True, False
                func = "exonic" if g[4] != "nonessential" else "ncRNA_intronic"
                aachange = "-"
            else:
                start = stop = pos_left = pos
                ref_c, alt_c, ref_aa, alt_aa = CODON_TABLE[klass]
                ref = ref_c[0].upper()
                var = alt_c[0].upper()
                is_snv, is_indel, is_mnv = True, False, False
                func = "exonic"
                aachange = f"{g[0]}:NM_XXX:exon1:c.{ref}1{var}:p.{ref_aa}1{alt_aa}"
            muts.append(
                dict(
                    clone_name=clone.name,
                    chr=g[1], start=start, stop=stop, pos_left=pos_left,
                    ref=ref, var=var,
                    is_snv=is_snv, is_indel=is_indel, is_mnv=is_mnv,
                    gene=g[0], gene_ensembl=g[3], func=func, aachange=aachange,
                    variant_class=klass,
                )
            )
    return muts


# Cache per-tumour mutation lists so muttable and AlphaMissense share them.
_MUTATIONS_CACHE: dict[str, list[dict]] = {}


def _mutations_for(spec: TumourSpec) -> list[dict]:
    key = spec.tumour_id_canonical
    if key not in _MUTATIONS_CACHE:
        _MUTATIONS_CACHE[key] = _assign_mutations(spec)
    return _MUTATIONS_CACHE[key]


def _write_muttable(specs: list[TumourSpec]) -> None:
    out = MOCK_ROOT / "tx842_mutation_table.tsv.gz"
    out.parent.mkdir(parents=True, exist_ok=True)

    # Full header from the catalog spec.
    header = [
        "patient_id", "patient_tumour", "chr", "start", "stop",
        "sample_name_hash", "ref", "var", "pos",
        "var_count", "ref_count", "cov", "vaf",
        "max_var_count", "min_depth", "max_vaf", "is_present",
        "gl_sample_name_hash", "is_snv", "is_indel", "is_mnv",
        "varscan", "mutect", "vaf02",
        "varscan_any", "mutect_any", "vaf02_any",
        "gene", "func_refgene", "aachange_refgene",
        "snp137", "cosmic90_coding", "cadd13gt10",
        "ljb23_sift", "ljb23_pp2hvar", "gnomad211_af",
        "CALLED", "PASS",
        "pan_cancer_driver", "lung_driver", "luad_driver", "lusc_driver",
        "oncogene", "tumor_suppressor", "cosmic_count",
        "is_deleterious", "driver_category", "is_driver_mut",
        "tumour_id", "ffpe", "ffpe_confidence", "mutation_id",
        "cn_a", "cn_b", "major_cn", "minor_cn", "mutcpn",
        "phylo_ccf", "obs_ccf", "clone_region_proportion",
        "mean_cluster_region_ccf", "cluster_region_clonality",
        "is_trunk", "timing", "mutation_cluster",
        "clean_cluster", "tree_cluster", "tree_cn_removed_cluster",
        "is_cn_fail",
        "ith_state", "ith_state_ambiguity1",
        "is_homogen", "is_artefact_sig_fail", "is_homogen_ambiguity1",
        "rna_var_count", "rna_vaf", "rna_cov",
        "is_expressed", "is_expressed_any",
    ]

    rows = []
    for spec in specs:
        muts = _mutations_for(spec)
        sample_ids = [
            f"{spec.patient_id}_SU_T{spec.tumour_ordinal}-R{i + 1}--{RNG.randbytes(6).hex()}"
            for i in range(spec.n_samples)
        ]
        gl = f"{spec.patient_id}_BS_GL--{RNG.randbytes(6).hex()}"

        for m in muts:
            cluster_idx = int(m["clone_name"].removeprefix("clone"))
            mutation_id = f"{spec.patient_id}_{spec.tumour_ordinal}:{m['chr']}:{m['pos_left']}:{m['ref']}:{m['var']}"
            is_trunk = m["clone_name"] == "clone1"
            for s_idx, sid in enumerate(sample_ids):
                # Sample-level presence: 90% present in the "primary" sample per clone,
                # 60% detected in others via rediscovery.
                is_present = RNG.random() < (0.9 if s_idx == 0 else 0.6)
                vaf = round(RNG.uniform(0.05, 0.45), 4) if is_present else 0.0
                cov = RNG.randint(60, 350)
                var_count = int(vaf * cov)
                ref_count = cov - var_count

                # Local CN — coarse but plausible: WGD-carrying clones get higher CN.
                # Use `mutation_cluster` to look up cumulative WGD on the fly.
                cn_boost = _cumulative_gds(spec, m["clone_name"])
                major_cn = 2.0 + cn_boost
                minor_cn = 1.0 + max(0, cn_boost - 1)
                mutcpn = round(RNG.uniform(1.0, major_cn + minor_cn), 3)

                is_driver = m["gene"] in {"TP53", "KRAS", "EGFR", "MYC"}

                row = [
                    spec.patient_id, spec.tumour_id_muttable, m["chr"], m["start"], m["stop"],
                    sid, m["ref"], m["var"], m["pos_left"],
                    var_count, ref_count, cov, vaf,
                    max(var_count, 40), min(cov, 120), min(0.5, vaf + 0.05), is_present,
                    gl, m["is_snv"], m["is_indel"], m["is_mnv"],
                    True, RNG.random() < 0.7, vaf > 0.02,
                    True, True, True,
                    m["gene"], m["func"], m["aachange"],
                    "", "", m["variant_class"] == "nonsense",
                    "", "", "",
                    True,  # CALLED — set to True; we'll test False elsewhere in future
                    True,  # PASS — always True (matches Tx842)
                    is_driver, is_driver, is_driver and m["gene"] != "MYC",
                    is_driver and m["gene"] == "TP53",
                    m["gene"] in {"KRAS", "EGFR", "MYC"},
                    m["gene"] == "TP53",
                    0,
                    m["variant_class"] in ("nonsense", "indel"),
                    "1" if is_driver else "",
                    is_driver,
                    spec.tumour_ordinal, False, "high_conf",
                    mutation_id,
                    major_cn - minor_cn, minor_cn, major_cn, minor_cn, mutcpn,
                    round(RNG.uniform(0.6, 1.0), 3), round(RNG.uniform(0.6, 1.0), 3),
                    0.0, round(RNG.uniform(0.7, 1.0), 2),
                    "clonal" if is_trunk else "subclonal",
                    is_trunk, "early",
                    float(cluster_idx),           # mutation_cluster (float)
                    True, True, False, False,
                    1, 1, is_trunk, False, is_trunk,
                    "", "", "", False, False,
                ]
                rows.append(row)

    with gzip.open(out, "wt", newline="\n") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(x) for x in r) + "\n")

    print(f"wrote {out}  ({len(rows)} rows)")


def _cumulative_gds(spec: TumourSpec, clone_name: str) -> int:
    """Sum of gds_at_clone along the parent chain to root (inclusive)."""
    idx = {c.name: c for c in spec.clones}
    total = 0
    node: str | None = clone_name
    while node and node != "diploid":
        total += idx[node].gds_at_clone
        node = idx[node].parent
    return total


def _fmt(v: object) -> str:
    if v is None:
        return "NA"
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


# --- WGD calls ---------------------------------------------------------------

def _write_wgd_calls(specs: list[TumourSpec]) -> None:
    WGD_DIR.mkdir(parents=True, exist_ok=True)
    out = WGD_DIR / "WGD_calls.tsv"

    header = [
        "tumour_id", "clone", "parent",
        "number_of_gds_at_clone", "total_gds_per_tumour",
        "quality", "status", "class", "sub_class",
    ]

    lines = ["\t".join(header)]
    for spec in specs:
        for c in spec.clones:
            lines.append("\t".join([
                spec.tumour_id_canonical,
                c.name,
                c.parent,
                str(c.gds_at_clone),
                str(spec.total_gds),
                spec.quality,
                spec.status,
                spec.class_ if spec.status == "resolved" else "",
                spec.sub_class or "",
            ]))

    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


# --- Clinical ----------------------------------------------------------------

def _write_clinical(specs: list[TumourSpec]) -> None:
    patients: dict[str, list[TumourSpec]] = {}
    for s in specs:
        patients.setdefault(s.patient_id, []).append(s)

    header = [
        "REGTrialNo", "Patient_ID", "age", "sex", "ethnicity",
        "pre_surgery_FEV1", "registration_ECOG",
        "smoking_status", "smoking_status_group",
        "Number_cigs_perday_use", "Number_years_smoking_use",
        "pack_years_calculated", "pack_years_truncated",
        "age_start_smoking_survey", "age_start_smoking_exsmoker", "age_stop_smoking_exsmoker",
        "alcohol_units", "alcohol_frequency", "alcohol_frequency_over7units",
        "alcohol_years", "alcohol_years_truncated",
        "Surgery_type", "adjuvant_treatment_given", "adjuvant_treatment_YN",
        "Centrally_reviewed_histology_per_patient",
        "Lesion1_site_central.reviewed", "histology1_central.reviewed",
        "histology1_group_central.reviewed",
        "pTStage_v8_lesion1_central.reviewed", "pNStage_lesion1_central.reviewed",
        "pTNMStage_v8_lesion1_central.reviewed", "SizePath_lesion1_central.reviewed",
        "PathPleuInv_lesion1_central.reviewed", "margin_status_lesion1_central.reviewed",
        "Lesion2_site_central.reviewed", "histology2_central.reviewed",
        "histology2_group_central.reviewed",
        "pTStage_v8_lesion2_central.reviewed", "pNStage_lesion2_central.reviewed",
        "pTNMStage_v8_lesion2_central.reviewed", "SizePath_lesion2_central.reviewed",
        "PathPleuInv_lesion2_central.reviewed", "margin_status_lesion2_central.reviewed",
        "note_histology", "note_pathology_report",
        "new_primary_cancer", "new_primary_flag",
        "second_lesion_detected_in_baseline", "CTC_notes",
        "Relapse_cat_new", "Relapse_site_text",
        "Surgical_bed_recurrence", "Ipsilateral_lung", "Mediastinum",
        "Contralateral_lung", "Pleural_nodules_effusion",
        "Liver", "Adrenal", "Bone", "Brain",
        "Extrathoracic_lymph_nodes", "Other", "Brain_only",
        "extracranial_extrathoracic_single",
        "PrevCancerAge", "Previous_cancer",
        "HomeSmoke", "WorkSmoke30Yrs", "WorkArsenic", "WorkAsbestos",
        "WorkBenzene", "WorkBisphenolA", "WorkCH2O", "WorkChemicals",
        "WorkChromHex", "WorkDiesel", "WorkDioxins", "WorkMetals",
        "WorkNatFibres", "WorkPAHs", "WorkPBDEs", "WorkRadon",
        "WorkSolvents", "WorkVinylCl",
        "PathTissueFind",
    ]

    def _blank_row(pid: str) -> list[str]:
        return ["NA"] * len(header)

    lines = ["\t".join(header)]
    for pid, tumours in patients.items():
        row = _blank_row(pid)
        row[0] = f"A_{pid}"
        row[1] = pid
        row[2] = str(RNG.randint(45, 82))
        row[3] = RNG.choice(["Male", "Female"])
        row[4] = "White- British"
        row[7] = RNG.choice(["Never", "Ex-Smoker", "Current"])
        row[8] = row[7]
        row[11] = str(RNG.randint(0, 60))
        row[12] = row[11]
        row[21] = "Lobectomy"
        row[22] = "No adjuvant"
        row[23] = "No"
        row[24] = tumours[0].histology_group

        # Lesion1 always populated.
        t1 = tumours[0]
        row[25] = "Right Upper Lobe"
        row[26] = "Invasive adenocarcinoma" if t1.histology_group == "LUAD" else "Squamous cell carcinoma"
        row[27] = t1.histology_group
        row[28] = RNG.choice(["1a", "1b", "2a", "2b", "3a"])
        row[29] = "0"

        # Lesion2 populated iff multi-tumour.
        if len(tumours) > 1:
            t2 = tumours[1]
            row[34] = "Left Lower Lobe"
            row[35] = "Squamous cell carcinoma"
            row[36] = t2.histology_group
            row[37] = "1b"

        lines.append("\t".join(row))

    out = MOCK_ROOT / "tx842_clinical.tsv"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


# --- AlphaMissense (per-tumour VEP output) -----------------------------------

def _write_alphamissense(specs: list[TumourSpec]) -> None:
    AM_ROOT.mkdir(parents=True, exist_ok=True)

    header = [
        "#Uploaded_variation", "Location", "Allele", "Gene", "Feature",
        "Feature_type", "Consequence", "cDNA_position", "CDS_position",
        "Protein_position", "Amino_acids", "Codons", "Existing_variation",
        "IMPACT", "DISTANCE", "STRAND", "FLAGS", "am_class", "am_pathogenicity",
    ]
    comments = [
        "## ENSEMBL VARIANT EFFECT PREDICTOR v115.2 (mock)",
        "## assembly version GRCh37.p13",
        "## AlphaMissense plugin — thresholds benign<0.34, ambiguous, pathogenic>0.564",
    ]

    for spec in specs:
        muts = _mutations_for(spec)
        tdir = AM_ROOT / spec.tumour_id_canonical
        tdir.mkdir(parents=True, exist_ok=True)
        out = tdir / f"{spec.tumour_id_canonical}_muttable_annotated.tsv"

        lines = comments + ["\t".join(header)]
        upload_id = f"{spec.patient_id}_{spec.tumour_ordinal}"
        for m in muts:
            loc = f"{m['chr']}:{m['pos_left']}"
            if m["variant_class"] == "missense":
                am_score = round(RNG.uniform(0.6, 0.99) if m["gene"] in {"TP53", "KRAS", "EGFR"}
                                 else RNG.uniform(0.02, 0.35), 4)
                am_class = ("pathogenic" if am_score > 0.564
                            else "ambiguous" if am_score >= 0.34
                            else "benign")
                cons = "missense_variant"
                aa = m["aachange"].split("p.")[-1]
                aa_ref, aa_alt = aa[0], aa[-1]
                impact = "MODERATE"
                lines.append("\t".join([
                    upload_id, loc, m["var"], m["gene_ensembl"], "ENST00000000001",
                    "Transcript", cons, "100", "100", "34",
                    f"{aa_ref}/{aa_alt}", "Ctt/Gtt", "-", impact, "-", "1", "-",
                    am_class, str(am_score),
                ]))
            elif m["variant_class"] == "synonymous":
                lines.append("\t".join([
                    upload_id, loc, m["var"], m["gene_ensembl"], "ENST00000000001",
                    "Transcript", "synonymous_variant", "100", "100", "34",
                    "L/L", "Ctt/Ctc", "-", "LOW", "-", "1", "-", "-", "-",
                ]))
            elif m["variant_class"] == "nonsense":
                lines.append("\t".join([
                    upload_id, loc, m["var"], m["gene_ensembl"], "ENST00000000001",
                    "Transcript", "stop_gained", "100", "100", "34",
                    "R/*", "Cga/Tga", "-", "HIGH", "-", "1", "-", "-", "-",
                ]))
            else:  # indel
                lines.append("\t".join([
                    upload_id, loc, m["var"], m["gene_ensembl"], "ENST00000000001",
                    "Transcript", "intron_variant", "-", "-", "-",
                    "-", "-", "-", "MODIFIER", "-", "-1", "-", "-", "-",
                ]))
        out.write_text("\n".join(lines) + "\n")
    print(f"wrote AM files for {len(specs)} tumours under {AM_ROOT}")


# --- Driver list ------------------------------------------------------------

def _write_driver_list() -> None:
    """Emit a synthetic TX842-shaped driver gene list.

    Categorises our fixed gene panel:
      - `TP53`, `KRAS`, `EGFR` → lung_mut_driver TRUE
      - `MYC` → mut_driver TRUE but not lung-specific
      - others → not drivers
    """
    out_dir = MOCK_ROOT / "drivers"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "tx842_driver_gene_list.tsv"

    header = [
        "Hugo_Symbol", "gene_id", "driver_gene", "driver_role", "CGC_TIER",
        "mut_driver", "lung_mut_driver", "mut_driver_source",
        "CN_driver", "lung_CN_driver", "CNA_type", "CN_driver_source",
        "chr_hg38", "start_hg38", "end_hg38", "cyto_hg38", "strand_hg38",
        "chr_hg19", "start_hg19", "end_hg19", "cyto_hg19",
        "SYNONYMS",
    ]

    lung_mut_drivers = {"TP53", "KRAS", "EGFR"}
    pan_mut_drivers = lung_mut_drivers | {"MYC"}
    lung_cn_drivers = {"MYC", "EGFR"}
    cn_drivers = lung_cn_drivers | {"KRAS"}

    lines = ["\t".join(header)]
    for sym, chrom, start, ensembl, category in GENES:
        is_lung_mut = sym in lung_mut_drivers
        is_mut = sym in pan_mut_drivers
        is_lung_cn = sym in lung_cn_drivers
        is_cn = sym in cn_drivers
        driver = is_mut or is_cn

        role = {
            "TP53": "tumor_suppressor",
            "KRAS": "oncogene",
            "EGFR": "oncogene",
            "MYC":  "oncogene",
        }.get(sym, "NA")

        cgc_tier = "1" if sym in {"TP53", "KRAS", "EGFR", "MYC"} else "NA"
        cna_type = "amp" if is_cn and sym != "TP53" else ("NA" if not is_cn else "del")

        end = start + 990
        lines.append("\t".join([
            sym, ensembl,
            "TRUE" if driver else "FALSE",
            role,
            cgc_tier,
            "TRUE" if is_mut else "FALSE",
            "TRUE" if is_lung_mut else "FALSE",
            "CGC_v104" if is_mut else "NA",
            "TRUE" if is_cn else "FALSE",
            "TRUE" if is_lung_cn else "FALSE",
            cna_type,
            "CGC_v104" if is_cn else "NA",
            f"chr{chrom}", str(start), str(end), f"{chrom}q00.0", "+",
            f"chr{chrom}", str(start), str(end), f"{chrom}q00.0",
            "NA",
        ]))

    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


# --- Entry ------------------------------------------------------------------

def main() -> None:
    MOCK_ROOT.mkdir(parents=True, exist_ok=True)
    specs = cohort()
    print(f"cohort: {len(specs)} tumours across {len({s.patient_id for s in specs})} patients")

    _write_muttable(specs)
    _write_wgd_calls(specs)
    _write_clinical(specs)
    _write_alphamissense(specs)
    _write_driver_list()
    print("done.")


if __name__ == "__main__":
    main()
