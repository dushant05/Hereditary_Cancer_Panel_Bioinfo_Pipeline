#!/usr/bin/env python3

"""
Targeted NGS QC / Coverage Analysis

Inputs:
    1. BED file containing target regions
    2. BAM prefix OR BAM file

Example:
    python targeted_ngs_qc.py \
        XT-HS_sorted_hg38.bed \
        Patient_T

The script expects:
    Patient_T.bqsr.bam
    Patient_T.bqsr.bam.bai

Outputs:
    Patient_T.coverageGene.csv
    Patient_T.coverageExon.csv
    Patient_T.coverageChrY.csv
    Patient_T.QCSummary.csv
    Patient_T.QCSummary.pdf

QC criteria:
    Median Coverage       >= 700X
    On-target Reads       >= 80%
    High-quality Bases    >= 90%
    High-quality Reads    >= 85%
    N Bases               <= 2%
    Average Read Length   >= 100 bp
    Uniformity            reported as percentage

Dependencies:
    pip install pysam pandas numpy matplotlib reportlab
"""

import argparse
import os
import sys
import statistics
from collections import defaultdict

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.enums import TA_CENTER
from reportlab.platypus import (
    SimpleDocTemplate,
    Table,
    TableStyle,
    Paragraph,
    Spacer,
    Image,
)
from reportlab.lib.units import inch


# ============================================================
# CONFIGURATION
# ============================================================

MEDIAN_COVERAGE_CUTOFF = 700
ONTARGET_CUTOFF = 80.0
HQ_BASE_CUTOFF = 90.0
HQ_READ_CUTOFF = 85.0
N_BASE_CUTOFF = 2.0
READ_LENGTH_CUTOFF = 100

# Uniformity definition:
# percentage of target bases having coverage >= 0.20 * median coverage
UNIFORMITY_FRACTION = 0.20

# Base quality threshold
BASE_QUALITY = 20

# Mapping quality threshold
MAP_QUALITY = 20


# ============================================================
# LOGGING
# ============================================================

def log(message):
    print(message, flush=True)


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


# ============================================================
# ARGUMENTS
# ============================================================

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Targeted NGS BAM/BED QC and coverage analysis"
    )

    parser.add_argument(
        "bed",
        help="Target BED file"
    )

    parser.add_argument(
        "bam_prefix",
        help="BAM prefix, e.g. Patient_T OR complete BAM filename"
    )

    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Output prefix. Default is BAM prefix without .bqsr.bam"
    )

    return parser.parse_args()


# ============================================================
# FILE HANDLING
# ============================================================

def resolve_bam(bam_prefix):

    if os.path.isfile(bam_prefix):
        return bam_prefix

    if os.path.isfile(bam_prefix + ".bqsr.bam"):
        return bam_prefix + ".bqsr.bam"

    fail(
        f"Cannot find BAM file.\n"
        f"Tried:\n"
        f"  {bam_prefix}\n"
        f"  {bam_prefix}.bqsr.bam"
    )


def resolve_bam_index(bam):

    possible_indexes = [
        bam + ".bai",
        bam.replace(".bam", ".bai")
    ]

    for index in possible_indexes:
        if os.path.isfile(index):
            return index

    fail(
        f"BAM index not found for {bam}.\n"
        f"Please run:\n"
        f"samtools index {bam}"
    )


# ============================================================
# BED READING
# ============================================================

def read_bed(bed_file):

    log(f"Reading BED file: {bed_file}")

    columns = [
        "chrom",
        "start",
        "end",
        "gene"
    ]

    records = []

    with open(bed_file) as fh:

        for line_number, line in enumerate(fh, 1):

            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                continue

            chrom = fields[0]

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue

            gene = fields[3] if len(fields) >= 4 else "Unknown"

            records.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "gene": gene
                }
            )

    if not records:
        fail("No valid regions found in BED file.")

    bed = pd.DataFrame(records)

    bed["length"] = bed["end"] - bed["start"]

    return bed


# ============================================================
# BAM STATISTICS
# ============================================================

def calculate_read_statistics(bam):

    log("Calculating read statistics...")

    total_reads = 0
    mapped_reads = 0
    duplicate_reads = 0
    secondary_reads = 0
    supplementary_reads = 0

    total_bases = 0
    high_quality_bases = 0
    n_bases = 0

    read_lengths = []

    with pysam.AlignmentFile(bam, "rb") as sam:

        for read in sam.fetch(until_eof=True):

            if read.is_unmapped:
                continue

            if read.is_secondary:
                secondary_reads += 1
                continue

            if read.is_supplementary:
                supplementary_reads += 1
                continue

            total_reads += 1
            mapped_reads += 1

            if read.is_duplicate:
                duplicate_reads += 1

            length = read.query_length

            if length:
                read_lengths.append(length)

            qualities = read.query_qualities

            if qualities is None:
                continue

            sequence = read.query_sequence

            total_bases += len(qualities)

            high_quality_bases += sum(
                q >= BASE_QUALITY for q in qualities
            )

            if sequence:
                n_bases += sum(
                    base.upper() == "N"
                    for base in sequence
                )

    if total_reads == 0:
        fail("No mapped reads found in BAM.")

    average_read_length = (
        np.mean(read_lengths)
        if read_lengths
        else 0
    )

    high_quality_base_percent = (
        100.0 * high_quality_bases / total_bases
        if total_bases
        else 0
    )

    n_base_percent = (
        100.0 * n_bases / total_bases
        if total_bases
        else 0
    )

    return {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "duplicate_reads": duplicate_reads,
        "secondary_reads": secondary_reads,
        "supplementary_reads": supplementary_reads,
        "total_bases": total_bases,
        "high_quality_bases": high_quality_bases,
        "high_quality_base_percent": high_quality_base_percent,
        "n_bases": n_bases,
        "n_base_percent": n_base_percent,
        "average_read_length": average_read_length
    }


# ============================================================
# ON-TARGET READS
# ============================================================

def calculate_ontarget_reads(bam, bed):

    log("Calculating on-target reads...")

    total_primary_mapped = 0
    on_target_reads = set()

    with pysam.AlignmentFile(bam, "rb") as sam:

        for read in sam.fetch(until_eof=True):

            if read.is_unmapped:
                continue

            if read.is_secondary:
                continue

            if read.is_supplementary:
                continue

            total_primary_mapped += 1

            if read.mapping_quality < MAP_QUALITY:
                continue

            try:
                blocks = read.get_blocks()
            except Exception:
                continue

            for chrom, start, end in bed[
                ["chrom", "start", "end"]
            ].itertuples(index=False):

                if read.reference_name != chrom:
                    continue

                # Quick overlap test
                overlap = False

                for block_start, block_end in blocks:

                    if (
                        block_start < end
                        and block_end > start
                    ):
                        overlap = True
                        break

                if overlap:
                    on_target_reads.add(read.query_name)
                    break

    on_target_percent = (
        100.0 * len(on_target_reads) / total_primary_mapped
        if total_primary_mapped
        else 0
    )

    return {
        "total_primary_mapped": total_primary_mapped,
        "on_target_reads": len(on_target_reads),
        "on_target_percent": on_target_percent
    }


# ============================================================
# TARGET COVERAGE
# ============================================================

def calculate_target_coverage(bam, bed):

    log("Calculating target coverage...")

    rows = []

    with pysam.AlignmentFile(bam, "rb") as sam:

        for index, region in bed.iterrows():

            chrom = region["chrom"]
            start = int(region["start"])
            end = int(region["end"])
            gene = region["gene"]

            try:
                coverage = np.array(
                    sam.count_coverage(
                        chrom,
                        start,
                        end,
                        quality_threshold=0
                    )
                ).sum(axis=0)

            except ValueError:
                coverage = np.zeros(end - start)

            if len(coverage) == 0:
                continue

            avg_cov = float(np.mean(coverage))
            median_cov = float(np.median(coverage))

            q20 = np.mean(coverage >= 20) * 100
            q100 = np.mean(coverage >= 100) * 100
            q200 = np.mean(coverage >= 200) * 100
            q500 = np.mean(coverage >= 500) * 100
            q700 = np.mean(coverage >= 700) * 100

            rows.append(
                {
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "gene": gene,
                    "length": end - start,
                    "avgCoverage": round(avg_cov, 2),
                    "medianCoverage": round(median_cov, 2),
                    "percent_20X": round(q20, 2),
                    "percent_100X": round(q100, 2),
                    "percent_200X": round(q200, 2),
                    "percent_500X": round(q500, 2),
                    "percent_700X": round(q700, 2)
                }
            )

    if not rows:
        fail("Unable to calculate coverage for any BED region.")

    return pd.DataFrame(rows)


# ============================================================
# GENE COVERAGE
# ============================================================

def calculate_gene_coverage(exon_coverage):

    log("Calculating gene-level coverage...")

    grouped = (
        exon_coverage
        .groupby("gene", dropna=False)
        .agg(
            avgCoverage=("avgCoverage", "mean"),
            medianCoverage=("medianCoverage", "median"),
            percent_20X=("percent_20X", "mean"),
            percent_100X=("percent_100X", "mean"),
            percent_200X=("percent_200X", "mean"),
            percent_500X=("percent_500X", "mean"),
            percent_700X=("percent_700X", "mean")
        )
        .reset_index()
    )

    numeric_columns = [
        "avgCoverage",
        "medianCoverage",
        "percent_20X",
        "percent_100X",
        "percent_200X",
        "percent_500X",
        "percent_700X"
    ]

    grouped[numeric_columns] = grouped[numeric_columns].round(2)

    return grouped


# ============================================================
# GLOBAL TARGET COVERAGE
# ============================================================

def calculate_global_coverage(bam, bed):

    log("Calculating global target coverage...")

    all_coverage = []

    with pysam.AlignmentFile(bam, "rb") as sam:

        for _, region in bed.iterrows():

            chrom = region["chrom"]
            start = int(region["start"])
            end = int(region["end"])

            try:

                cov = np.array(
                    sam.count_coverage(
                        chrom,
                        start,
                        end,
                        quality_threshold=0
                    )
                ).sum(axis=0)

                all_coverage.extend(cov.tolist())

            except ValueError:
                continue

    if not all_coverage:
        fail("No target coverage could be calculated.")

    cov = np.array(all_coverage)

    mean_cov = float(np.mean(cov))
    median_cov = float(np.median(cov))

    uniformity_threshold = median_cov * UNIFORMITY_FRACTION

    uniformity = (
        np.mean(cov >= uniformity_threshold) * 100
    )

    coverage_700 = np.mean(cov >= 700) * 100
    coverage_500 = np.mean(cov >= 500) * 100
    coverage_200 = np.mean(cov >= 200) * 100
    coverage_100 = np.mean(cov >= 100) * 100
    coverage_20 = np.mean(cov >= 20) * 100

    return {
        "mean_coverage": mean_cov,
        "median_coverage": median_cov,
        "uniformity_threshold": uniformity_threshold,
        "uniformity_percent": uniformity,
        "percent_20X": coverage_20,
        "percent_100X": coverage_100,
        "percent_200X": coverage_200,
        "percent_500X": coverage_500,
        "percent_700X": coverage_700,
        "coverage_array": cov
    }


# ============================================================
# CHR Y COVERAGE
# ============================================================

def calculate_chrY_coverage(exon_coverage):

    chr_y = exon_coverage[
        exon_coverage["chromosome"].isin(
            ["chrY", "Y"]
        )
    ].copy()

    return chr_y


# ============================================================
# QC PASS / FAIL
# ============================================================

def qc_status(value, cutoff, direction):

    if direction == ">=":
        return "PASS" if value >= cutoff else "FAIL"

    if direction == "<=":
        return "PASS" if value <= cutoff else "FAIL"

    return "NA"


def generate_qc_summary(
    global_cov,
    read_stats,
    on_target
):

    median_cov = global_cov["median_coverage"]

    ontarget = on_target["on_target_percent"]

    hq_bases = read_stats["high_quality_base_percent"]

    # High-quality reads:
    # Reads having >= 85% bases with Q >= 20
    # This is calculated separately below in the BAM pass.
    hq_reads = calculate_high_quality_reads

    summary = [
        {
            "Parameter": "Median Coverage",
            "Result": round(median_cov, 2),
            "Criterion": ">= 700X",
            "Remark": qc_status(
                median_cov,
                MEDIAN_COVERAGE_CUTOFF,
                ">="
            )
        },
        {
            "Parameter": "On-target Reads",
            "Result": round(ontarget, 2),
            "Criterion": ">= 80%",
            "Remark": qc_status(
                ontarget,
                ONTARGET_CUTOFF,
                ">="
            )
        },
        {
            "Parameter": "High-quality Bases",
            "Result": round(hq_bases, 2),
            "Criterion": ">= 90%",
            "Remark": qc_status(
                hq_bases,
                HQ_BASE_CUTOFF,
                ">="
            )
        },
        {
            "Parameter": "High-quality Reads",
            "Result": round(hq_reads, 2),
            "Criterion": ">= 85%",
            "Remark": qc_status(
                hq_reads,
                HQ_READ_CUTOFF,
                ">="
            )
        },
        {
            "Parameter": "N Bases",
            "Result": round(
                read_stats["n_base_percent"],
                2
            ),
            "Criterion": "<= 2%",
            "Remark": qc_status(
                read_stats["n_base_percent"],
                N_BASE_CUTOFF,
                "<="
            )
        },
        {
            "Parameter": "Average Read Length",
            "Result": round(
                read_stats["average_read_length"],
                2
            ),
            "Criterion": ">= 100 bp",
            "Remark": qc_status(
                read_stats["average_read_length"],
                READ_LENGTH_CUTOFF,
                ">="
            )
        },
        {
            "Parameter": "Uniformity",
            "Result": round(
                global_cov["uniformity_percent"],
                2
            ),
            "Criterion": ">= 80% of target bases >= 20% median",
            "Remark": "REPORT"
        }
    ]

    return pd.DataFrame(summary)


# ============================================================
# HIGH QUALITY READS
# ============================================================

def calculate_high_quality_reads(bam):

    log("Calculating high-quality reads...")

    total_reads = 0
    high_quality_reads = 0

    with pysam.AlignmentFile(bam, "rb") as sam:

        for read in sam.fetch(until_eof=True):

            if read.is_unmapped:
                continue

            if read.is_secondary:
                continue

            if read.is_supplementary:
                continue

            if read.query_qualities is None:
                continue

            total_reads += 1

            qualities = np.array(
                read.query_qualities
            )

            if len(qualities) == 0:
                continue

            fraction_q20 = np.mean(
                qualities >= BASE_QUALITY
            ) * 100

            if fraction_q20 >= 85:
                high_quality_reads += 1

    if total_reads == 0:
        return 0

    return (
        high_quality_reads / total_reads
    ) * 100


# ============================================================
# PDF HELPERS
# ============================================================

def make_table(dataframe, font_size=7):

    data = [
        list(dataframe.columns)
    ]

    for row in dataframe.itertuples(index=False):
        data.append(
            [
                str(x)
                for x in row
            ]
        )

    table = Table(
        data,
        repeatRows=1,
        hAlign="LEFT"
    )

    table.setStyle(
        TableStyle(
            [
                (
                    "BACKGROUND",
                    (0, 0),
                    (-1, 0),
                    colors.HexColor("#2F5597")
                ),
                (
                    "TEXTCOLOR",
                    (0, 0),
                    (-1, 0),
                    colors.white
                ),
                (
                    "FONTNAME",
                    (0, 0),
                    (-1, 0),
                    "Helvetica-Bold"
                ),
                (
                    "FONTSIZE",
                    (0, 0),
                    (-1, -1),
                    font_size
                ),
                (
                    "GRID",
                    (0, 0),
                    (-1, -1),
                    0.25,
                    colors.grey
                ),
                (
                    "VALIGN",
                    (0, 0),
                    (-1, -1),
                    "MIDDLE"
                ),
                (
                    "ROWBACKGROUNDS",
                    (0, 1),
                    (-1, -1),
                    [
                        colors.white,
                        colors.HexColor("#F2F2F2")
                    ]
                )
            ]
        )
    )

    return table


# ============================================================
# COVERAGE HISTOGRAM
# ============================================================

def create_coverage_plot(
    coverage_array,
    output_file
):

    plt.figure(
        figsize=(8, 5)
    )

    # Limit plot to avoid extreme high-coverage values
    plot_values = coverage_array[
        coverage_array <= np.percentile(
            coverage_array,
            99
        )
    ]

    plt.hist(
        plot_values,
        bins=50,
        color="lightblue",
        edgecolor="black"
    )

    plt.xlabel("Coverage")
    plt.ylabel("Number of target bases")
    plt.title("Target Coverage Distribution")

    plt.axvline(
        MEDIAN_COVERAGE_CUTOFF,
        color="red",
        linestyle="--",
        label="700X QC threshold"
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=150
    )

    plt.close()


# ============================================================
# PDF REPORT
# ============================================================

def create_pdf_report(
    output_pdf,
    sample,
    qc_summary,
    read_stats,
    on_target,
    global_cov,
    coverage_plot
):

    log(f"Creating QC PDF: {output_pdf}")

    styles = getSampleStyleSheet()

    title_style = styles["Title"]
    title_style.alignment = TA_CENTER

    doc = SimpleDocTemplate(
        output_pdf,
        pagesize=letter,
        rightMargin=30,
        leftMargin=30,
        topMargin=30,
        bottomMargin=30
    )

    story = []

    story.append(
        Paragraph(
            f"Targeted NGS QC Summary<br/>{sample}",
            title_style
        )
    )

    story.append(
        Spacer(1, 15)
    )

    # --------------------------------------------------------
    # QC summary
    # --------------------------------------------------------

    story.append(
        Paragraph(
            "QC Summary",
            styles["Heading2"]
        )
    )

    story.append(
        make_table(
            qc_summary,
            font_size=8
        )
    )

    story.append(
        Spacer(1, 15)
    )

    # --------------------------------------------------------
    # Read statistics
    # --------------------------------------------------------

    read_df = pd.DataFrame(
        [
            ["Total mapped reads",
             read_stats["mapped_reads"]],

            ["Duplicate reads",
             read_stats["duplicate_reads"]],

            ["Total bases",
             read_stats["total_bases"]],

            ["High-quality bases (%)",
             round(
                 read_stats[
                     "high_quality_base_percent"
                 ],
                 2
             )],

            ["N bases (%)",
             round(
                 read_stats["n_base_percent"],
                 2
             )],

            ["Average read length (bp)",
             round(
                 read_stats[
                     "average_read_length"
                 ],
                 2
             )],

            ["On-target reads",
             on_target["on_target_reads"]],

            ["On-target reads (%)",
             round(
                 on_target[
                     "on_target_percent"
                 ],
                 2
             )]
        ],
        columns=[
            "Parameter",
            "Value"
        ]
    )

    story.append(
        Paragraph(
            "Read Statistics",
            styles["Heading2"]
        )
    )

    story.append(
        make_table(
            read_df,
            font_size=8
        )
    )

    story.append(
        Spacer(1, 15)
    )

    # --------------------------------------------------------
    # Coverage statistics
    # --------------------------------------------------------

    coverage_df = pd.DataFrame(
        [
            ["Mean coverage",
             round(
                 global_cov[
                     "mean_coverage"
                 ],
                 2
             )],

            ["Median coverage",
             round(
                 global_cov[
                     "median_coverage"
                 ],
                 2
             )],

            ["20X or greater (%)",
             round(
                 global_cov[
                     "percent_20X"
                 ],
                 2
             )],

            ["100X or greater (%)",
             round(
                 global_cov[
                     "percent_100X"
                 ],
                 2
             )],

            ["200X or greater (%)",
             round(
                 global_cov[
                     "percent_200X"
                 ],
                 2
             )],

            ["500X or greater (%)",
             round(
                 global_cov[
                     "percent_500X"
                 ],
                 2
             )],

            ["700X or greater (%)",
             round(
                 global_cov[
                     "percent_700X"
                 ],
                 2
             )],

            [
                "Uniformity (%)",
                round(
                    global_cov[
                        "uniformity_percent"
                    ],
                    2
                )
            ],

            [
                "Uniformity threshold",
                f"{global_cov['uniformity_threshold']:.2f}X"
            ]
        ],
        columns=[
            "Coverage Metric",
            "Value"
        ]
    )

    story.append(
        Paragraph(
            "Coverage Statistics",
            styles["Heading2"]
        )
    )

    story.append(
        make_table(
            coverage_df,
            font_size=8
        )
    )

    story.append(
        Spacer(1, 15)
    )

    # --------------------------------------------------------
    # Histogram
    # --------------------------------------------------------

    story.append(
        Paragraph(
            "Coverage Distribution",
            styles["Heading2"]
        )
    )

    story.append(
        Image(
            coverage_plot,
            width=7 * inch,
            height=4.2 * inch
        )
    )

    doc.build(story)


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_arguments()

    bed_file = os.path.abspath(args.bed)

    bam_file = resolve_bam(
        args.bam_prefix
    )

    resolve_bam_index(
        bam_file
    )

    if args.output_prefix:

        output_prefix = args.output_prefix

    else:

        output_prefix = os.path.basename(
            bam_file
        )

        if output_prefix.endswith(
            ".bqsr.bam"
        ):
            output_prefix = output_prefix[
                :-len(".bqsr.bam")
            ]

        elif output_prefix.endswith(
            ".bam"
        ):
            output_prefix = output_prefix[
                :-len(".bam")
            ]

    log("")
    log("=" * 70)
    log("TARGETED NGS QC ANALYSIS")
    log("=" * 70)
    log(f"BED : {bed_file}")
    log(f"BAM : {bam_file}")
    log(f"OUT : {output_prefix}")
    log("=" * 70)

    # --------------------------------------------------------
    # BED
    # --------------------------------------------------------

    bed = read_bed(
        bed_file
    )

    # --------------------------------------------------------
    # Read statistics
    # --------------------------------------------------------

    read_stats = calculate_read_statistics(
        bam_file
    )

    # --------------------------------------------------------
    # High quality reads
    # --------------------------------------------------------

    global calculate_high_quality_reads
    calculate_high_quality_reads = calculate_high_quality_reads(
        bam_file
    )

    hq_reads_percent = calculate_high_quality_reads

    # --------------------------------------------------------
    # On target
    # --------------------------------------------------------

    on_target = calculate_ontarget_reads(
        bam_file,
        bed
    )

    # --------------------------------------------------------
    # Exon/target coverage
    # --------------------------------------------------------

    exon_coverage = calculate_target_coverage(
        bam_file,
        bed
    )

    # --------------------------------------------------------
    # Gene coverage
    # --------------------------------------------------------

    gene_coverage = calculate_gene_coverage(
        exon_coverage
    )

    # --------------------------------------------------------
    # Global coverage
    # --------------------------------------------------------

    global_cov = calculate_global_coverage(
        bam_file,
        bed
    )

    # --------------------------------------------------------
    # chrY
    # --------------------------------------------------------

    chr_y = calculate_chrY_coverage(
        exon_coverage
    )

    # --------------------------------------------------------
    # QC summary
    # --------------------------------------------------------

    qc_summary = pd.DataFrame(
        [
            {
                "Parameter": "Median Coverage",
                "Result": round(
                    global_cov["median_coverage"],
                    2
                ),
                "Criterion": ">= 700X",
                "Remark": qc_status(
                    global_cov["median_coverage"],
                    MEDIAN_COVERAGE_CUTOFF,
                    ">="
                )
            },
            {
                "Parameter": "On-target Reads",
                "Result": round(
                    on_target["on_target_percent"],
                    2
                ),
                "Criterion": ">= 80%",
                "Remark": qc_status(
                    on_target["on_target_percent"],
                    ONTARGET_CUTOFF,
                    ">="
                )
            },
            {
                "Parameter": "High-quality Bases",
                "Result": round(
                    read_stats[
                        "high_quality_base_percent"
                    ],
                    2
                ),
                "Criterion": ">= 90%",
                "Remark": qc_status(
                    read_stats[
                        "high_quality_base_percent"
                    ],
                    HQ_BASE_CUTOFF,
                    ">="
                )
            },
            {
                "Parameter": "High-quality Reads",
                "Result": round(
                    hq_reads_percent,
                    2
                ),
                "Criterion": ">= 85%",
                "Remark": qc_status(
                    hq_reads_percent,
                    HQ_READ_CUTOFF,
                    ">="
                )
            },
            {
                "Parameter": "N Bases",
                "Result": round(
                    read_stats[
                        "n_base_percent"
                    ],
                    2
                ),
                "Criterion": "<= 2%",
                "Remark": qc_status(
                    read_stats[
                        "n_base_percent"
                    ],
                    N_BASE_CUTOFF,
                    "<="
                )
            },
            {
                "Parameter": "Average Read Length",
                "Result": round(
                    read_stats[
                        "average_read_length"
                    ],
                    2
                ),
                "Criterion": ">= 100 bp",
                "Remark": qc_status(
                    read_stats[
                        "average_read_length"
                    ],
                    READ_LENGTH_CUTOFF,
                    ">="
                )
            },
            {
                "Parameter": "Uniformity",
                "Result": round(
                    global_cov[
                        "uniformity_percent"
                    ],
                    2
                ),
                "Criterion": ">= 80%",
                "Remark": (
                    "PASS"
                    if global_cov[
                        "uniformity_percent"
                    ] >= 80
                    else "FAIL"
                )
            }
        ]
    )

    # --------------------------------------------------------
    # Write CSV outputs
    # --------------------------------------------------------

    log("Writing CSV files...")

    gene_coverage.to_csv(
        f"{output_prefix}.coverageGene.csv",
        index=False
    )

    exon_coverage.to_csv(
        f"{output_prefix}.coverageExon.csv",
        index=False
    )

    chr_y.to_csv(
        f"{output_prefix}.coverageChrY.csv",
        index=False
    )

    qc_summary.to_csv(
        f"{output_prefix}.QCSummary.csv",
        index=False
    )

    # --------------------------------------------------------
    # Coverage plot
    # --------------------------------------------------------

    coverage_plot = (
        f"{output_prefix}.coverage_histogram.png"
    )

    create_coverage_plot(
        global_cov["coverage_array"],
        coverage_plot
    )

    # --------------------------------------------------------
    # PDF
    # --------------------------------------------------------

    output_pdf = (
        f"{output_prefix}.QCSummary.pdf"
    )

    create_pdf_report(
        output_pdf,
        output_prefix,
        qc_summary,
        read_stats,
        on_target,
        global_cov,
        coverage_plot
    )

    # --------------------------------------------------------
    # Final console report
    # --------------------------------------------------------

    log("")
    log("=" * 70)
    log("QC SUMMARY")
    log("=" * 70)

    for _, row in qc_summary.iterrows():

        log(
            f"{row['Parameter']:<25} "
            f"{row['Result']:<12} "
            f"{row['Criterion']:<30} "
            f"{row['Remark']}"
        )

    log("=" * 70)

    log("")
    log("Output files:")

    log(
        f"  {output_prefix}.coverageGene.csv"
    )

    log(
        f"  {output_prefix}.coverageExon.csv"
    )

    log(
        f"  {output_prefix}.coverageChrY.csv"
    )

    log(
        f"  {output_prefix}.QCSummary.csv"
    )

    log(
        f"  {output_prefix}.QCSummary.pdf"
    )

    log("")
    log("QC analysis completed successfully.")


if __name__ == "__main__":
    main()
