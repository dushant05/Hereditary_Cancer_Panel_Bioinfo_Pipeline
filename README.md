# Hereditary Cancer Panel Pipeline

A reproducible bioinformatics pipeline for the analysis of hereditary cancer gene panel sequencing data, from raw sequencing reads through variant identification, annotation, and clinical interpretation. The pipeline is designed to streamline the processing of targeted NGS data and provide a standardized workflow for detecting germline single-nucleotide variants (SNVs), small insertions/deletions (indels), and other clinically relevant variants across genes associated with hereditary cancer predisposition.

## Key Features
1. End-to-end processing of targeted NGS data
2. Raw FASTQ quality control and preprocessing
3. Read alignment to the reference genome
4. BAM processing, sorting, and quality assessment
5. Germline variant calling
6. Variant quality filtering and prioritization
7. Variant annotation using relevant population and clinical databases
8. Identification of potentially pathogenic and clinically relevant variants
9. Generation of analysis-ready variant reports
10. Reproducible and modular workflow suitable for research and diagnostic environments

## Workflow

Raw FASTQ
   │
   ▼
Quality Control
   │
   ▼
Read Trimming / Preprocessing
   │
   ▼
Reference Genome Alignment
   │
   ▼
BAM Processing & QC
   │
   ▼
Germline Variant Calling
   │
   ▼
Variant Filtering
   │
   ▼
Variant Annotation
   │
   ▼
Variant Prioritization
   │
   ▼
Clinical / Research Report


## Requirements

The analysis server must have:

- Bash
- fastp
- FastQC
- BWA
- SAMtools
- Picard
- GATK
- R
- VEP

Reference data:

- hg38 reference genome
- Hereditary Cancer panel BED
- Gene BED
- dbSNP
- Mills and 1000G indel database
- VEP databases

## Quick Start
For a new server:

cd hereditary-cancer-panel

cp config/config.sh.example config/config.sh

nano config/config.sh

chmod +x scripts/run.sh

chmod +x scripts/hereditary_cancer.sh

# Run the patient:

./scripts/run.sh \
    PATIENT001 \
    /data/fastq/PATIENT001_R1.fastq.gz \
    /data/fastq/PATIENT001_R2.fastq.gz

# For the coverage QC check, run:

python3 coverage_qc.py \
    --bam  \
    --bed  \
    --output 

Monitor:

tail -f /data/hereditary-cancer/results/PATIENT001/pipeline.log

## Pipeline Output
The main outputs are:

1. Filtered VCF
2. VEP annotated VCF
3. CNV CSV
4. QC reports
5. BAM
6. BAM index
7. GATK recalibration files
8. Pipeline logs

## Intended Use

This pipeline is intended to support the analysis of hereditary cancer panel sequencing data for research, method development, and bioinformatics workflows. Final clinical interpretation and reporting should be performed according to applicable laboratory procedures, validated protocols, and relevant clinical guidelines.

## Reproducibility
The workflow is organized into modular steps to facilitate reproducible analysis, parameter tracking, and integration into automated NGS processing environments. Configuration files and workflow parameters can be adapted for different sequencing platforms, gene panels, reference genomes, and annotation resources.

## Pipeline Version
Hereditary Cancer Panel Pipeline
Version: 1.0.0
Reference: GRCh38 / hg38
