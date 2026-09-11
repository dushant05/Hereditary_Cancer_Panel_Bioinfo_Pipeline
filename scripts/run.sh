#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# Hereditary Cancer Panel
# Pipeline Launcher
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

PIPELINE="${SCRIPT_DIR}/hereditary_cancer.sh"
CONFIG="${PROJECT_DIR}/config/config.sh"

echo "============================================================"
echo "        HEREDITARY CANCER PANEL ANALYSIS"
echo "============================================================"

# ------------------------------------------------------------
# Check configuration
# ------------------------------------------------------------

if [[ ! -f "$CONFIG" ]]; then
    echo ""
    echo "ERROR: config.sh not found."
    echo ""
    echo "Create it using:"
    echo ""
    echo "cp config/config.sh.example config/config.sh"
    echo ""
    exit 1
fi

source "$CONFIG"

# ------------------------------------------------------------
# Check pipeline
# ------------------------------------------------------------

if [[ ! -f "$PIPELINE" ]]; then
    echo "ERROR: Pipeline not found:"
    echo "$PIPELINE"
    exit 1
fi

chmod +x "$PIPELINE"

# ------------------------------------------------------------
# Input arguments
# ------------------------------------------------------------

if [[ $# -ne 3 ]]; then

    echo ""
    echo "Usage:"
    echo ""
    echo "  ./scripts/run.sh <PATIENT_ID> <R1_FASTQ> <R2_FASTQ>"
    echo ""
    echo "Example:"
    echo ""
    echo "  ./scripts/run.sh \\"
    echo "      PATIENT001 \\"
    echo "      /data/fastq/PATIENT001_R1.fastq.gz \\"
    echo "      /data/fastq/PATIENT001_R2.fastq.gz"
    echo ""

    exit 1
fi

PATIENT_ID="$1"
R1="$2"
R2="$3"

# ------------------------------------------------------------
# Validate input files
# ------------------------------------------------------------

if [[ ! -f "$R1" ]]; then
    echo "ERROR: R1 FASTQ does not exist:"
    echo "$R1"
    exit 1
fi

if [[ ! -f "$R2" ]]; then
    echo "ERROR: R2 FASTQ does not exist:"
    echo "$R2"
    exit 1
fi

# ------------------------------------------------------------
# Create result directory
# ------------------------------------------------------------

RESULT_DIR="${RESULTS_DIR}/${PATIENT_ID}"

mkdir -p "$RESULT_DIR"

LOG="${RESULT_DIR}/launcher.log"

# ------------------------------------------------------------
# Start logging
# ------------------------------------------------------------

echo "==================================================" | tee -a "$LOG"
echo "Hereditary Cancer Panel" | tee -a "$LOG"
echo "==================================================" | tee -a "$LOG"

echo "Patient ID : $PATIENT_ID" | tee -a "$LOG"
echo "R1         : $R1" | tee -a "$LOG"
echo "R2         : $R2" | tee -a "$LOG"
echo "Threads    : $THREADS" | tee -a "$LOG"
echo "Start      : $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$LOG"

# ------------------------------------------------------------
# Run pipeline
# ------------------------------------------------------------

echo "" | tee -a "$LOG"
echo "Starting pipeline..." | tee -a "$LOG"

"$PIPELINE" \
    "$PATIENT_ID" \
    "$R1" \
    "$R2" \
    "$RESULT_DIR" \
    2>&1 | tee "${RESULT_DIR}/pipeline.log"

STATUS=${PIPESTATUS[0]}

# ------------------------------------------------------------
# Finish
# ------------------------------------------------------------

if [[ "$STATUS" -eq 0 ]]; then

    echo "" | tee -a "$LOG"
    echo "==================================================" | tee -a "$LOG"
    echo "ANALYSIS COMPLETED SUCCESSFULLY" | tee -a "$LOG"
    echo "Patient : $PATIENT_ID" | tee -a "$LOG"
    echo "Results : $RESULT_DIR" | tee -a "$LOG"
    echo "End     : $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$LOG"
    echo "==================================================" | tee -a "$LOG"

else

    echo "" | tee -a "$LOG"
    echo "==================================================" | tee -a "$LOG"
    echo "ANALYSIS FAILED" | tee -a "$LOG"
    echo "Patient : $PATIENT_ID" | tee -a "$LOG"
    echo "Check   : ${RESULT_DIR}/pipeline.log" | tee -a "$LOG"
    echo "==================================================" | tee -a "$LOG"

    exit "$STATUS"

fi

