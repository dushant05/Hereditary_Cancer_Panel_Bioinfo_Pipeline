# Hereditary Cancer Panel Pipeline

NGS analysis pipeline for hereditary cancer panel sequencing.

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
