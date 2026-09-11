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

## Pipeline

FASTQ
→ FASTP
→ FastQC
→ BWA MEM
→ SAMtools
→ Picard
→ MarkDuplicates
→ GATK BQSR
→ GATK HaplotypeCaller
→ SNP/INDEL filtering
→ CNV analysis
→ VEP annotation
→ Final results

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

## Configuration

Copy the example configuration:

```bash
cp config/config.sh.example config/config.sh

# Hereditary_Cancer_Panel_Bioinfo_Pipeline
