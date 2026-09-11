#!/usr/bin/env bash

# ============================================================
# Hereditary Cancer Panel Configuration
# ============================================================

# CPU
THREADS=24


# ============================================================
# REFERENCE GENOME
# ============================================================

REFERENCE="/data/reference/hg38/hg38.fasta"


# ============================================================
# PANEL
# ============================================================

PANEL_BED="/data/panel/Hereditary_Cancer_panel.bed"

GENES_BED="/data/panel/GCP_Gene_hg38.bed"


# ============================================================
# GATK DATABASES
# ============================================================

MILLS_INDELS="/data/database/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

DBSNP="/data/database/dbsnp_146.hg38.vcf.gz"


# ============================================================
# CNV
# ============================================================

CNV_SCRIPT="/opt/hereditary-cancer/GCP_call_CNV_single_sample.R"


# ============================================================
# VEP
# ============================================================

VEP_PATH="/opt/ensembl-vep"

VEP_DATA="/data/vep"


# ============================================================
# READ GROUP
# ============================================================

RG_LIBRARY="library"

RG_PLATFORM="illumina"

RG_PLATFORM_UNIT="HaloPlex"


# ============================================================
# RESULTS
# ============================================================

RESULTS_DIR="/data/hereditary-cancer/results"
