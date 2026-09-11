#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Hereditary Cancer Panel
# Main Analysis Pipeline
# ============================================================

if [[ $# -ne 4 ]]; then
    echo "Usage:"
    echo "  $0 <PATIENT_ID> <R1_FASTQ> <R2_FASTQ> <RESULT_DIR>"
    exit 1
fi

patientID="$1"
rawFqTR1="$2"
rawFqTR2="$3"
RESULT_DIR="$4"

mkdir -p "$RESULT_DIR"

cd "$RESULT_DIR"

echo "==================================================" >> output.log
echo "HEREDITARY CANCER PANEL" >> output.log
echo "==================================================" >> output.log

echo "PATIENT: ${patientID}" >> output.log
echo "START: $(date '+%Y-%m-%d %H:%M:%S')" >> output.log

echo "R1: ${rawFqTR1}" >> output.log
echo "R2: ${rawFqTR2}" >> output.log

# ============================================================
# CONFIGURATION
# ============================================================

echo "Loading pipeline configuration..." >> output.log

# These variables come from config.sh
#
# REFERENCE
# PANEL_BED
# GENES_BED
# MILLS_INDELS
# DBSNP
# CNV_SCRIPT
# VEP_PATH
# VEP_DATA
# THREADS

# ============================================================
# DIRECTORIES
# ============================================================

mkdir -p clean
mkdir -p QC
mkdir -p Output

# ============================================================
# YOUR ORIGINAL PIPELINE STARTS HERE
# ============================================================

echo "ANALYSIS STARTED" >> output.log
date +"%Y-%m-%d %H:%M:%S" >> output.log


# ------------------------------------------------------------
# FASTP
# ------------------------------------------------------------

echo "Starting FASTP..." >> output.log

fastp \
    -h "clean/${patientID}_T.fastp.html" \
    -j "clean/${patientID}_T.fastp.json" \
    -w "$THREADS" \
    -l 30 \
    --cut_mean_quality 20 \
    -a auto \
    -f 10 \
    -F 10 \
    -p \
    -z 9 \
    -i "$rawFqTR1" \
    -I "$rawFqTR2" \
    -o "clean/${patientID}_T_1.clean.fastq.gz" \
    -O "clean/${patientID}_T_2.clean.fastq.gz"


# ------------------------------------------------------------
# FASTQC
# ------------------------------------------------------------

echo "Starting FASTQC..." >> output.log

fastqc \
    -t 2 \
    "clean/${patientID}_T_1.clean.fastq.gz" \
    "clean/${patientID}_T_2.clean.fastq.gz" \
    -o QC/


# ------------------------------------------------------------
# ALIGNMENT
# ------------------------------------------------------------

echo "Starting BWA alignment..." >> output.log

bwa mem \
    -M \
    -t "$THREADS" \
    "$REFERENCE" \
    "clean/${patientID}_T_1.clean.fastq.gz" \
    "clean/${patientID}_T_2.clean.fastq.gz" |
    samtools view -bh |
    samtools sort \
        -@ "$THREADS" \
        -o "${patientID}_T.alnpe.sort.bam"

samtools index "${patientID}_T.alnpe.sort.bam"


# ------------------------------------------------------------
# ADD READ GROUP
# ------------------------------------------------------------

echo "Adding read groups..." >> output.log

picard AddOrReplaceReadGroups \
    I="${patientID}_T.alnpe.sort.bam" \
    O="${patientID}_T.sortrg.bam" \
    RGLB="$RG_LIBRARY" \
    RGPL="$RG_PLATFORM" \
    RGPU="$RG_PLATFORM_UNIT" \
    RGSM="$patientID" \
    COMPRESSION_LEVEL=0 \
    VALIDATION_STRINGENCY=LENIENT \
    SO=coordinate \
    CREATE_INDEX=true


# ------------------------------------------------------------
# MARK DUPLICATES
# ------------------------------------------------------------

echo "Marking duplicates..." >> output.log

picard MarkDuplicates \
    I="${patientID}_T.sortrg.bam" \
    O="${patientID}_T.sort.markdup.bam" \
    REMOVE_DUPLICATES=true \
    AS=true \
    METRICS_FILE="${patientID}_T.sort.markdup.metrics"

samtools index "${patientID}_T.sort.markdup.bam"


# ------------------------------------------------------------
# BQSR
# ------------------------------------------------------------

echo "Running BQSR..." >> output.log

gatk BaseRecalibrator \
    -I "${patientID}_T.sort.markdup.bam" \
    -R "$REFERENCE" \
    --known-sites "$MILLS_INDELS" \
    --known-sites "$DBSNP" \
    -O "${patientID}_T.recal_data.table"


gatk ApplyBQSR \
    -R "$REFERENCE" \
    -I "${patientID}_T.sort.markdup.bam" \
    --bqsr-recal-file "${patientID}_T.recal_data.table" \
    -O "${patientID}_T.bqsr.bam"

samtools index "${patientID}_T.bqsr.bam"


# ------------------------------------------------------------
# HAPLOTYPECALLER
# ------------------------------------------------------------

echo "Running GATK HaplotypeCaller..." >> output.log

gatk HaplotypeCaller \
    -L "$PANEL_BED" \
    -R "$REFERENCE" \
    -I "${patientID}_T.bqsr.bam" \
    -O "${patientID}.haplotypeCaller_raw.vcf" \
    -stand-call-conf 30 \
    --dbsnp "$DBSNP"


# ------------------------------------------------------------
# SNP
# ------------------------------------------------------------

gatk SelectVariants \
    -R "$REFERENCE" \
    -V "${patientID}.haplotypeCaller_raw.vcf" \
    --select-type-to-include SNP \
    -O "${patientID}.haplotypeCaller_raw_snps.vcf"


gatk VariantFiltration \
    -R "$REFERENCE" \
    -V "${patientID}.haplotypeCaller_raw_snps.vcf" \
    --filter-expression "(QD < 2.0) || (FS > 60.0) || (MQ < 40.0) || (MQRankSum < -12.5) || (ReadPosRankSum < -8.0) || (SOR > 4.0)" \
    --filter-name "${patientID}.snp.filter" \
    -O "${patientID}.haplotypeCaller_filtered_snps.vcf"


# ------------------------------------------------------------
# INDEL
# ------------------------------------------------------------

gatk SelectVariants \
    -R "$REFERENCE" \
    -V "${patientID}.haplotypeCaller_raw.vcf" \
    --select-type-to-include INDEL \
    -O "${patientID}.haplotypeCaller_raw_indels.vcf"


gatk VariantFiltration \
    -R "$REFERENCE" \
    -V "${patientID}.haplotypeCaller_raw_indels.vcf" \
    --filter-expression "(QD < 2.0) || (FS > 200.0) || (ReadPosRankSum < -20.0) || (SOR > 10.0)" \
    --filter-name "${patientID}.indel.filter" \
    -O "${patientID}.haplotypeCaller_filtered_indels.vcf"


# ------------------------------------------------------------
# COMBINE
# ------------------------------------------------------------

gatk SortVcf \
    -I "${patientID}.haplotypeCaller_filtered_snps.vcf" \
    -I "${patientID}.haplotypeCaller_filtered_indels.vcf" \
    -O "${patientID}.haplotypeCaller_filter_marked.vcf"


grep -i '^#\|PASS' \
    "${patientID}.haplotypeCaller_filter_marked.vcf" \
    > "${patientID}.haplotypeCaller_filtered.vcf"


# ------------------------------------------------------------
# CNV
# ------------------------------------------------------------

echo "Running CNV analysis..." >> output.log

Rscript "$CNV_SCRIPT" \
    "$GENES_BED" \
    "${patientID}_T"


# ------------------------------------------------------------
# VEP
# ------------------------------------------------------------

echo "Running VEP..." >> output.log

export VEP_PATH
export VEP_DATA
export PERL5LIB="${VEP_PATH}:${PERL5LIB:-}"


"$VEP_PATH/vep" \
    --species homo_sapiens \
    --assembly GRCh38 \
    --offline \
    --no_stats \
    --sift b \
    --ccds \
    --uniprot \
    --hgvs \
    --symbol \
    --numbers \
    --domains \
    --gene_phenotype \
    --canonical \
    --protein \
    --biotype \
    --tsl \
    --pubmed \
    --variant_class \
    --shift_hgvs 1 \
    --check_existing \
    --total_length \
    --allele_number \
    --no_escape \
    --xref_refseq \
    --failed 1 \
    --vcf \
    --minimal \
    --flag_pick_allele \
    --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir "$VEP_DATA" \
    --fasta "$REFERENCE" \
    --input_file "${patientID}.haplotypeCaller_filtered.vcf" \
    --output_file "${patientID}_T.haplotypeCaller.vep.vcf" \
    --polyphen b \
    --af \
    --af_1kg \
    --regulatory


# ------------------------------------------------------------
# FINISHED
# ------------------------------------------------------------

echo "==================================================" >> output.log
echo "DATA PROCESSING COMPLETED" >> output.log
echo "PATIENT: ${patientID}" >> output.log
echo "END: $(date '+%Y-%m-%d %H:%M:%S')" >> output.log
echo "==================================================" >> output.log

echo ""
echo "Analysis completed for ${patientID}"
echo "Results: ${RESULT_DIR}"
