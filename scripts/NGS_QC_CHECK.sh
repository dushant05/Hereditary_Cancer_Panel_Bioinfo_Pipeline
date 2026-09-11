#!/usr/bin/env bash

# ============================================================
# NGS Alignment QC
# Illumina paired-end targeted sequencing
#
# Usage:
#   ./qc.sh sample.bam targets.bed sample_R1.fastq.gz sample_R2.fastq.gz
#
# FASTQ files are optional:
#   ./qc.sh sample.bam targets.bed
#
# Requirements:
#   samtools
#   bedtools
#   awk
#   gzip
#
# Output:
#   sample_QC.tsv
#   sample_flagstat.txt
#   sample_target_depth.txt
#   sample_target_coverage.txt
# ============================================================

set -euo pipefail

# -----------------------------
# Check arguments
# -----------------------------
if [[ $# -lt 2 ]]; then
    echo ""
    echo "Usage:"
    echo "  $0 sample.bam targets.bed [R1.fastq.gz] [R2.fastq.gz]"
    echo ""
    exit 1
fi

BAM="$1"
BED="$2"

R1="${3:-}"
R2="${4:-}"

# -----------------------------
# Check input files
# -----------------------------
if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM file not found: $BAM"
    exit 1
fi

if [[ ! -f "$BED" ]]; then
    echo "ERROR: BED file not found: $BED"
    exit 1
fi

# -----------------------------
# Check dependencies
# -----------------------------
for tool in samtools bedtools awk; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: $tool is not installed or not in PATH"
        exit 1
    fi
done

# -----------------------------
# Sample name
# -----------------------------
SAMPLE=$(basename "$BAM")
SAMPLE="${SAMPLE%.bam}"
SAMPLE="${SAMPLE%.cram}"

echo ""
echo "=========================================="
echo "        NGS ALIGNMENT QC"
echo "=========================================="
echo "Sample : $SAMPLE"
echo "BAM    : $BAM"
echo "BED    : $BED"
echo "=========================================="
echo ""

# -----------------------------
# Output files
# -----------------------------
QC="${SAMPLE}_QC.tsv"
FLAGSTAT="${SAMPLE}_flagstat.txt"
DEPTH="${SAMPLE}_target_depth.txt"
COVERAGE="${SAMPLE}_target_coverage.txt"

# -----------------------------
# BAM index check
# -----------------------------
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
    echo "BAM index not found."
    echo "Creating BAM index..."
    samtools index "$BAM"
fi

# -----------------------------
# 1. samtools flagstat
# -----------------------------
echo "[1/7] Calculating alignment statistics..."

samtools flagstat "$BAM" > "$FLAGSTAT"

TOTAL_READS=$(awk '/in total/ {print $1; exit}' "$FLAGSTAT")

MAPPED_READS=$(awk '/ mapped \\(/ {
    gsub("\\(", "", $5);
    print $1;
    exit
}' "$FLAGSTAT")

PROPERLY_PAIRED=$(awk '/properly paired/ {
    print $1;
    exit
}' "$FLAGSTAT")

DUPLICATES=$(awk '/duplicates/ {
    print $1;
    exit
}' "$FLAGSTAT")

# Avoid empty values
TOTAL_READS=${TOTAL_READS:-0}
MAPPED_READS=${MAPPED_READS:-0}
PROPERLY_PAIRED=${PROPERLY_PAIRED:-0}
DUPLICATES=${DUPLICATES:-0}

# Mapping %
if [[ "$TOTAL_READS" -gt 0 ]]; then
    MAPPING_RATE=$(awk -v m="$MAPPED_READS" -v t="$TOTAL_READS" \
        'BEGIN {printf "%.2f",100*m/t}')
else
    MAPPING_RATE="0.00"
fi

# Properly paired %
if [[ "$TOTAL_READS" -gt 0 ]]; then
    PROPER_PAIR_RATE=$(awk -v p="$PROPERLY_PAIRED" -v t="$TOTAL_READS" \
        'BEGIN {printf "%.2f",100*p/t}')
else
    PROPER_PAIR_RATE="0.00"
fi

# Duplicate %
if [[ "$MAPPED_READS" -gt 0 ]]; then
    DUP_RATE=$(awk -v d="$DUPLICATES" -v m="$MAPPED_READS" \
        'BEGIN {printf "%.2f",100*d/m}')
else
    DUP_RATE="0.00"
fi

# -----------------------------
# 2. Target BED cleanup
# -----------------------------
echo "[2/7] Preparing target BED..."

TARGET_BED="${SAMPLE}_targets.cleaned.bed"

awk 'BEGIN{OFS="\t"}
     $1 !~ /^#/ && NF >= 3 {
         if ($2 >= 0 && $3 > $2)
             print $1,$2,$3
     }' "$BED" > "$TARGET_BED"

# -----------------------------
# 3. On-target reads
# -----------------------------
echo "[3/7] Calculating on-target reads..."

# Count aligned reads overlapping target regions.
# -f 0.1 requires at least 10% of each read to overlap a target.
ON_TARGET_READS=$(
    samtools view -c -F 4 -L "$TARGET_BED" "$BAM"
)

if [[ "$MAPPED_READS" -gt 0 ]]; then
    ON_TARGET_RATE=$(awk \
        -v o="$ON_TARGET_READS" \
        -v m="$MAPPED_READS" \
        'BEGIN {printf "%.2f",100*o/m}')
else
    ON_TARGET_RATE="0.00"
fi

# -----------------------------
# 4. Target depth
# -----------------------------
echo "[4/7] Calculating target depth..."

# samtools depth:
# -a = include zero coverage positions
#
# Restrict to target BED.
samtools depth -a -b "$TARGET_BED" "$BAM" > "$DEPTH"

# Mean depth
MEAN_DEPTH=$(
    awk '
    {
        sum += $3
        n++
    }
    END {
        if (n > 0)
            printf "%.2f", sum/n
        else
            print "0.00"
    }' "$DEPTH"
)

# Median depth
MEDIAN_DEPTH=$(
    awk '
    {
        depth[++n] = $3
    }
    END {
        if (n == 0) {
            print "0.00"
            exit
        }

        # sort
        for (i=1; i<=n; i++) {
            for (j=i+1; j<=n; j++) {
                if (depth[i] > depth[j]) {
                    tmp=depth[i]
                    depth[i]=depth[j]
                    depth[j]=tmp
                }
            }
        }

        if (n % 2 == 1)
            printf "%.2f", depth[(n+1)/2]
        else
            printf "%.2f", (depth[n/2]+depth[n/2+1])/2
    }' "$DEPTH"
)

# -----------------------------
# 5. Coverage thresholds
# -----------------------------
echo "[5/7] Calculating coverage thresholds..."

TOTAL_TARGET_BASES=$(wc -l < "$DEPTH")

if [[ "$TOTAL_TARGET_BASES" -eq 0 ]]; then
    echo "ERROR: No target bases found."
    exit 1
fi

calc_coverage() {
    local threshold="$1"

    awk -v t="$threshold" '
    {
        if ($3 >= t)
            n++
    }
    END {
        if (NR > 0)
            printf "%.2f",100*n/NR
        else
            print "0.00"
    }' "$DEPTH"
}

COV_1X=$(calc_coverage 1)
COV_10X=$(calc_coverage 10)
COV_20X=$(calc_coverage 20)
COV_30X=$(calc_coverage 30)
COV_50X=$(calc_coverage 50)
COV_100X=$(calc_coverage 100)
COV_500X=$(calc_coverage 500)
COV_1000X=$(calc_coverage 1000)

# -----------------------------
# 6. Uniformity
# -----------------------------
echo "[6/7] Calculating uniformity..."

# Uniformity definition used here:
#
#   % target bases with depth between
#   0.2 x mean depth and 2 x mean depth
#
# This is one commonly used definition.
#
LOW=$(awk -v m="$MEAN_DEPTH" 'BEGIN {printf "%.6f",0.2*m}')
HIGH=$(awk -v m="$MEAN_DEPTH" 'BEGIN {printf "%.6f",2*m}')

UNIFORMITY=$(
    awk -v low="$LOW" -v high="$HIGH" '
    {
        if ($3 >= low && $3 <= high)
            n++
    }
    END {
        if (NR > 0)
            printf "%.2f",100*n/NR
        else
            print "0.00"
    }' "$DEPTH"
)

# Alternative common metric:
# % bases >= 0.2 x mean depth
UNIFORMITY_20P=$(
    awk -v low="$LOW" '
    {
        if ($3 >= low)
            n++
    }
    END {
        if (NR > 0)
            printf "%.2f",100*n/NR
        else
            print "0.00"
    }' "$DEPTH"
)

# -----------------------------
# 7. Q30 from FASTQ
# -----------------------------
echo "[7/7] Calculating Q30..."

Q30="NA"

if [[ -n "$R1" && -n "$R2" ]]; then

    if [[ ! -f "$R1" ]]; then
        echo "WARNING: R1 FASTQ not found: $R1"
    elif [[ ! -f "$R2" ]]; then
        echo "WARNING: R2 FASTQ not found: $R2"
    else

        # Calculate fraction of bases with Phred >= 30.
        #
        # FASTQ quality encoding:
        # ASCII character - 33 = Phred score
        #
        # awk reads every 4th line (quality line).
        #
        Q30=$(
            {
                zcat "$R1"
                zcat "$R2"
            } | awk '
            NR % 4 == 0 {
                for (i=1; i<=length($0); i++) {
                    q = index("!\"#$%&'\''()*+,-./0123456789:;<=>?@ABCDEFGHIJ",
                              substr($0,i,1)) - 1

                    total++

                    if (q >= 30)
                        q30++
                }
            }

            END {
                if (total > 0)
                    printf "%.2f",100*q30/total
                else
                    print "0.00"
            }'
        )
    fi
fi

# -----------------------------
# Coverage summary by target
# -----------------------------
bedtools coverage \
    -a "$TARGET_BED" \
    -b "$BAM" \
    -mean > "$COVERAGE"

# -----------------------------
# Write final QC table
# -----------------------------
echo ""
echo "Writing QC report..."

cat > "$QC" <<EOF
Sample	$SAMPLE
Total_Reads	$TOTAL_READS
Mapped_Reads	$MAPPED_READS
Mapping_Rate_%	$MAPPING_RATE
Properly_Paired_Reads	$PROPERLY_PAIRED
Properly_Paired_%	$PROPER_PAIR_RATE
Duplicate_Reads	$DUPLICATES
Duplicate_Rate_%	$DUP_RATE
On_Target_Reads	$ON_TARGET_READS
On_Target_Rate_%	$ON_TARGET_RATE
Mean_Target_Depth	$MEAN_DEPTH
Median_Target_Depth	$MEDIAN_DEPTH
Target_Bases	$TOTAL_TARGET_BASES
Coverage_1X_%	$COV_1X
Coverage_10X_%	$COV_10X
Coverage_20X_%	$COV_20X
Coverage_30X_%	$COV_30X
Coverage_50X_%	$COV_50X
Coverage_100X_%	$COV_100X
Coverage_500X_%	$COV_500X
Coverage_1000X_%	$COV_1000X
Uniformity_0.2x-2xMean_%	$UNIFORMITY
Uniformity_>=0.2xMean_%	$UNIFORMITY_20P
Q30_%	$Q30
EOF

# -----------------------------
# Display results
# -----------------------------
echo ""
echo "=========================================="
echo "             QC SUMMARY"
echo "=========================================="

column -t -s $'\t' "$QC" 2>/dev/null || cat "$QC"

echo ""
echo "=========================================="
echo "Output files:"
echo "  $QC"
echo "  $FLAGSTAT"
echo "  $DEPTH"
echo "  $COVERAGE"
echo "=========================================="
echo ""

# Cleanup temporary BED
rm -f "$TARGET_BED"
