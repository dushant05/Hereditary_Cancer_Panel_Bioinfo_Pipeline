Hereditary Cancer Panel NGS Pipeline

Automated NGS analysis pipeline for a hereditary cancer panel.

Pipeline Workflow
FASTQ
  │
  ▼
FASTP
  │
  ▼
FastQC
  │
  ▼
BWA MEM
  │
  ▼
SAMtools
  │
  ▼
Picard AddOrReplaceReadGroups
  │
  ▼
Picard MarkDuplicates
  │
  ▼
GATK BaseRecalibrator
  │
  ▼
GATK ApplyBQSR
  │
  ▼
GATK HaplotypeCaller
  │
  ├── SNP
  │
  └── INDEL
       │
       ▼
   Variant Filtering
       │
       ▼
      CNV
       │
       ▼
  CNVPanelizer
       │
       ▼
      VEP
       │
       ▼
FINAL OUTPUT

1. Repository Structure
hereditary-cancer-panel/
│
├── .github/
│   └── workflows/
│       └── run-pipeline.yml
│
├── config/
│   ├── config.sh.example
│   └── config.sh
│
├── scripts/
│   ├── run.sh
│   ├── hereditary_cancer.sh
│   └── GCP_call_CNV_single_sample.R
│
├── README.md
├── .gitignore
└── LICENSE

2. Software Requirements

The analysis server should have:

Bash
fastp
FastQC
BWA
SAMtools
Picard
GATK
R
CNVPanelizer
VEP
Perl
GNU Parallel


Check the programs:

bash --version
fastp --version
fastqc --version
bwa
samtools --version
picard -h
gatk --version
R --version
perl --version
parallel --version


Check VEP:

vep --help


If VEP is installed in a custom directory:

/opt/ensembl-vep/vep --help

3. Clone the Repository

Clone the GitHub repository:

git clone YOUR_GITHUB_REPOSITORY


Enter the repository:

cd hereditary-cancer-panel


Check the files:

ls -lah


You should see:

README.md
scripts/
config/
.github/
.gitignore
LICENSE

4. Configure the Pipeline

Copy the example configuration:

cp config/config.sh.example config/config.sh


Open the configuration:

nano config/config.sh


Example:

#!/usr/bin/env bash

THREADS=24

REFERENCE="/data/reference/hg38/hg38.fasta"

PANEL_BED="/data/panel/Hereditary_Cancer_panel.bed"

GENES_BED="/data/panel/GCP_Gene_hg38.bed"

MILLS_INDELS="/data/database/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

DBSNP="/data/database/dbsnp_146.hg38.vcf.gz"

CNV_SCRIPT="/opt/hereditary-cancer/scripts/GCP_call_CNV_single_sample.R"

CNV_NORMAL_DIR="/data/Normal_data"

VEP_PATH="/opt/ensembl-vep"

VEP_DATA="/data/vep"

RG_LIBRARY="library"

RG_PLATFORM="illumina"

RG_PLATFORM_UNIT="HaloPlex"

RESULTS_DIR="/data/hereditary-cancer/results"


Save the file.

5. Check Reference Files

Check the reference genome:

ls -lh /data/reference/hg38/hg38.fasta


Check the FASTA index:

ls -lh /data/reference/hg38/hg38.fasta.fai


Check the GATK dictionary:

ls -lh /data/reference/hg38/hg38.dict


Check the panel BED:

ls -lh /data/panel/Hereditary_Cancer_panel.bed


Check the CNV gene BED:

ls -lh /data/panel/GCP_Gene_hg38.bed


Check dbSNP:

ls -lh /data/database/dbsnp_146.hg38.vcf.gz


Check dbSNP index:

ls -lh /data/database/dbsnp_146.hg38.vcf.gz.tbi


Check Mills:

ls -lh /data/database/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz


Check Mills index:

ls -lh /data/database/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

6. Check CNV Normal Samples

The CNV pipeline requires normal reference BAM files.

Check the directory:

ls -lh /data/Normal_data/


Check only BQSR BAM files:

find /data/Normal_data/ -type f -name "*bqsr.bam"


Count normal samples:

find /data/Normal_data/ -type f -name "*bqsr.bam" | wc -l


Check BAM indexes:

find /data/Normal_data/ -type f -name "*bqsr.bam.bai"


Example:

/data/Normal_data/Normal01.bqsr.bam
/data/Normal_data/Normal01.bqsr.bam.bai
/data/Normal_data/Normal02.bqsr.bam
/data/Normal_data/Normal02.bqsr.bam.bai
/data/Normal_data/Normal03.bqsr.bam
/data/Normal_data/Normal03.bqsr.bam.bai

7. Install/Check CNVPanelizer

Start R:

R


Inside R:

library(CNVPanelizer)


If it loads successfully:

packageVersion("CNVPanelizer")


Exit R:

q()

8. Make Scripts Executable

Run:

chmod +x scripts/run.sh


Run:

chmod +x scripts/hereditary_cancer.sh


Run:

chmod +x scripts/GCP_call_CNV_single_sample.R


Check:

ls -lh scripts/


You should see executable permissions such as:

-rwxr-xr-x run.sh
-rwxr-xr-x hereditary_cancer.sh
-rwxr-xr-x GCP_call_CNV_single_sample.R

9. Check the Pipeline Before Running

Run the launcher without arguments:

./scripts/run.sh


You should get:

Usage:

./scripts/run.sh <PATIENT_ID> <R1_FASTQ> <R2_FASTQ>

10. Run One Patient

Example:

./scripts/run.sh \
    PATIENT001 \
    /data/fastq/PATIENT001_R1.fastq.gz \
    /data/fastq/PATIENT001_R2.fastq.gz


Or in one line:

./scripts/run.sh PATIENT001 /data/fastq/PATIENT001_R1.fastq.gz /data/fastq/PATIENT001_R2.fastq.gz

11. Example Input

Input files:

/data/fastq/PATIENT001_R1.fastq.gz
/data/fastq/PATIENT001_R2.fastq.gz


Run:

./scripts/run.sh \
    PATIENT001 \
    /data/fastq/PATIENT001_R1.fastq.gz \
    /data/fastq/PATIENT001_R2.fastq.gz

12. Analysis Steps

The pipeline performs:

FASTP
fastp \
    -i PATIENT001_R1.fastq.gz \
    -I PATIENT001_R2.fastq.gz \
    -o PATIENT001_T_1.clean.fastq.gz \
    -O PATIENT001_T_2.clean.fastq.gz

FastQC
fastqc \
    PATIENT001_T_1.clean.fastq.gz \
    PATIENT001_T_2.clean.fastq.gz

BWA Alignment
bwa mem \
    -M \
    -t 24 \
    /data/reference/hg38/hg38.fasta \
    PATIENT001_T_1.clean.fastq.gz \
    PATIENT001_T_2.clean.fastq.gz

SAMtools Sorting
samtools sort \
    -@ 24 \
    -o PATIENT001_T.alnpe.sort.bam

BAM Index
samtools index PATIENT001_T.alnpe.sort.bam

Picard Read Groups
picard AddOrReplaceReadGroups \
    I=PATIENT001_T.alnpe.sort.bam \
    O=PATIENT001_T.sortrg.bam \
    RGLB=library \
    RGPL=illumina \
    RGPU=HaloPlex \
    RGSM=PATIENT001 \
    SO=coordinate \
    CREATE_INDEX=true

Picard MarkDuplicates
picard MarkDuplicates \
    I=PATIENT001_T.sortrg.bam \
    O=PATIENT001_T.sort.markdup.bam \
    REMOVE_DUPLICATES=true \
    AS=true \
    METRICS_FILE=PATIENT001_T.sort.markdup.metrics

GATK BQSR
gatk BaseRecalibrator \
    -I PATIENT001_T.sort.markdup.bam \
    -R /data/reference/hg38/hg38.fasta \
    --known-sites /data/database/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --known-sites /data/database/dbsnp_146.hg38.vcf.gz \
    -O PATIENT001_T.recal_data.table


Then:

gatk ApplyBQSR \
    -R /data/reference/hg38/hg38.fasta \
    -I PATIENT001_T.sort.markdup.bam \
    --bqsr-recal-file PATIENT001_T.recal_data.table \
    -O PATIENT001_T.bqsr.bam

GATK HaplotypeCaller
gatk HaplotypeCaller \
    -L /data/panel/Hereditary_Cancer_panel.bed \
    -R /data/reference/hg38/hg38.fasta \
    -I PATIENT001_T.bqsr.bam \
    -O PATIENT001.haplotypeCaller_raw.vcf \
    -stand-call-conf 30 \
    --dbsnp /data/database/dbsnp_146.hg38.vcf.gz

SNP Calling
gatk SelectVariants \
    -R /data/reference/hg38/hg38.fasta \
    -V PATIENT001.haplotypeCaller_raw.vcf \
    --select-type-to-include SNP \
    -O PATIENT001.haplotypeCaller_raw_snps.vcf

INDEL Calling
gatk SelectVariants \
    -R /data/reference/hg38/hg38.fasta \
    -V PATIENT001.haplotypeCaller_raw.vcf \
    --select-type-to-include INDEL \
    -O PATIENT001.haplotypeCaller_raw_indels.vcf

CNV

CNV is automatically called by:

Rscript \
    scripts/GCP_call_CNV_single_sample.R \
    /data/panel/GCP_Gene_hg38.bed \
    PATIENT001_T \
    /data/Normal_data


The CNV output is:

PATIENT001.copyNumberVariations.csv

VEP

VEP annotation is automatically performed on the filtered VCF.

13. Results

Results are stored in:

/data/hereditary-cancer/results/PATIENT001/


Example:

PATIENT001/
│
├── launcher.log
├── pipeline.log
├── output.log
│
├── clean/
│   ├── PATIENT001_T_1.clean.fastq.gz
│   ├── PATIENT001_T_2.clean.fastq.gz
│   ├── PATIENT001_T.fastp.html
│   └── PATIENT001_T.fastp.json
│
├── QC/
│   ├── FastQC reports
│   └── QC metrics
│
└── Output/
    ├── PATIENT001.haplotypeCaller_filtered.vcf
    ├── PATIENT001_T.haplotypeCaller.vep.vcf
    ├── PATIENT001.copyNumberVariations.csv
    └── other analysis files

14. View the Log

During analysis:

tail -f /data/hereditary-cancer/results/PATIENT001/pipeline.log


Or:

tail -f /data/hereditary-cancer/results/PATIENT001/output.log

15. Check the Final VCF
grep -v "^#" \
    /data/hereditary-cancer/results/PATIENT001/Output/PATIENT001.haplotypeCaller_filtered.vcf


Count variants:

grep -vc "^#" \
    /data/hereditary-cancer/results/PATIENT001/Output/PATIENT001.haplotypeCaller_filtered.vcf

16. Check CNV Results
cat \
    /data/hereditary-cancer/results/PATIENT001/Output/PATIENT001.copyNumberVariations.csv


Or:

column -s, -t \
    /data/hereditary-cancer/results/PATIENT001/Output/PATIENT001.copyNumberVariations.csv

17. Run Another Patient
./scripts/run.sh \
    PATIENT002 \
    /data/fastq/PATIENT002_R1.fastq.gz \
    /data/fastq/PATIENT002_R2.fastq.gz


Another:

./scripts/run.sh \
    PATIENT003 \
    /data/fastq/PATIENT003_R1.fastq.gz \
    /data/fastq/PATIENT003_R2.fastq.gz

18. Run Multiple Patients

Example:

./scripts/run.sh PATIENT001 /data/fastq/PATIENT001_R1.fastq.gz /data/fastq/PATIENT001_R2.fastq.gz
./scripts/run.sh PATIENT002 /data/fastq/PATIENT002_R1.fastq.gz /data/fastq/PATIENT002_R2.fastq.gz
./scripts/run.sh PATIENT003 /data/fastq/PATIENT003_R1.fastq.gz /data/fastq/PATIENT003_R2.fastq.gz


For production use, process samples according to the available CPU/RAM resources rather than launching many 24-thread jobs simultaneously.

19. Stop the Pipeline

If the pipeline needs to be stopped:

Ctrl+C

20. Update the Pipeline

Check GitHub for updates:

git pull


Check the current version:

git log --oneline -5


Check repository status:

git status

21. Upload Pipeline Changes to GitHub

After modifying the code:

git status


Add changes:

git add scripts/


Commit:

git commit -m "Update hereditary cancer pipeline"


Push:

git push

22. Important Data-Security Rule

Do NOT upload patient genomic data to GitHub.

Do NOT commit:

*.fastq
*.fastq.gz
*.bam
*.bai
*.cram
*.vcf
*.vcf.gz
*.csv


Do NOT commit:

config/config.sh


Only commit:

config/config.sh.example


Patient data should remain on the secure analysis server.

23. Quick Start

For a new server:

git clone YOUR_GITHUB_REPOSITORY

cd hereditary-cancer-panel

cp config/config.sh.example config/config.sh

nano config/config.sh

chmod +x scripts/run.sh

chmod +x scripts/hereditary_cancer.sh


Run the patient:

./scripts/run.sh \
    PATIENT001 \
    /data/fastq/PATIENT001_R1.fastq.gz \
    /data/fastq/PATIENT001_R2.fastq.gz


Monitor:

tail -f /data/hereditary-cancer/results/PATIENT001/pipeline.log

24. Pipeline Output

The main outputs are:

Filtered VCF
VEP annotated VCF
CNV CSV
QC reports
BAM
BAM index
GATK recalibration files
Pipeline logs

25. Pipeline Version
Hereditary Cancer Panel Pipeline
Version: 1.0.0
Reference: GRCh38 / hg38


Note: This pipeline is intended for research/analysis use unless it has been appropriately validated and authorized for clinical diagnostic use. Variant/CNV filtering parameters, reference datasets, software versions, and QC thresholds should be validated for the specific assay and laboratory.
