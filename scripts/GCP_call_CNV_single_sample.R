#!/usr/bin/env Rscript

############################################################
# Hereditary Cancer Panel
# CNV Calling - CNVPanelizer
############################################################

suppressPackageStartupMessages({
    library(CNVPanelizer)
})

############################################################
# ARGUMENTS
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    cat("\n")
    cat("Usage:\n")
    cat("  GCP_call_CNV_single_sample.R <GENES_BED> <SAMPLE_PREFIX> <NORMAL_DIR>\n")
    cat("\n")
    cat("Example:\n")
    cat("  GCP_call_CNV_single_sample.R \\\n")
    cat("      /data/panel/GCP_Gene_hg38.bed \\\n")
    cat("      PATIENT001_T \\\n")
    cat("      /data/Normal_data\n")
    cat("\n")
    quit(status = 1)
}

genesBed <- args[1]
samplePrefix <- args[2]
referenceDirectory <- args[3]

############################################################
# CHECK INPUTS
############################################################

if (!file.exists(genesBed)) {
    stop(
        paste(
            "ERROR: BED file does not exist:",
            genesBed
        )
    )
}

if (!dir.exists(referenceDirectory)) {
    stop(
        paste(
            "ERROR: Normal reference directory does not exist:",
            referenceDirectory
        )
    )
}

sampleBam <- paste0(samplePrefix, ".bqsr.bam")

if (!file.exists(sampleBam)) {
    stop(
        paste(
            "ERROR: Sample BAM does not exist:",
            sampleBam
        )
    )
}

############################################################
# START
############################################################

cat("\n")
cat("====================================================\n")
cat("HEREDITARY CANCER PANEL - CNV ANALYSIS\n")
cat("====================================================\n")

cat("BED file       :", genesBed, "\n")
cat("Sample prefix  :", samplePrefix, "\n")
cat("Sample BAM     :", sampleBam, "\n")
cat("Normal data    :", referenceDirectory, "\n")
cat("====================================================\n\n")


############################################################
# READING TARGETS
############################################################

cat("Reading target regions...\n")

genomicRangesFromBed <- BedToGenomicRanges(
    genesBed,
    ampliconColumn = 4,
    split = "_"
)

metadataFromGenomicRanges <- elementMetadata(
    genomicRangesFromBed
)

geneNames <- metadataFromGenomicRanges["geneNames"][, 1]

ampliconNames <- metadataFromGenomicRanges["ampliconNames"][, 1]


############################################################
# NORMAL REFERENCE SAMPLES
############################################################

cat("Reading normal reference BAM files...\n")

referenceFilenames <- list.files(
    path = referenceDirectory,
    pattern = "bqsr\\.bam$",
    full.names = TRUE
)

if (length(referenceFilenames) == 0) {
    stop(
        paste(
            "ERROR: No normal reference BAM files found in:",
            referenceDirectory
        )
    )
}

cat(
    "Number of normal reference samples:",
    length(referenceFilenames),
    "\n"
)

removePcrDuplicates <- FALSE

referenceReadCounts <- ReadCountsFromBam(
    referenceFilenames,
    genomicRangesFromBed,
    sampleNames = referenceFilenames,
    ampliconNames = ampliconNames,
    removeDup = removePcrDuplicates
)


############################################################
# SAMPLE READ COUNTS
############################################################

cat("Reading patient BAM file...\n")

sampleReadCounts <- ReadCountsFromBam(
    sampleBam,
    gr = genomicRangesFromBed,
    ampliconNames = ampliconNames,
    sampleNames = "Tumor",
    removeDup = FALSE
)


############################################################
# NORMALIZATION
############################################################

cat("Normalizing read counts...\n")

normalizedReadCounts <- CombinedNormalizedCounts(
    sampleReadCounts,
    referenceReadCounts,
    ampliconNames = ampliconNames
)


############################################################
# SPLIT NORMALIZED DATA
############################################################

samplesNormalizedReadCounts <-
    normalizedReadCounts["samples"][[1]]

referenceNormalizedReadCounts <-
    normalizedReadCounts["reference"][[1]]

cat(
    "Sample matrix dimensions:",
    dim(samplesNormalizedReadCounts),
    "\n"
)

cat(
    "Reference matrix dimensions:",
    dim(referenceNormalizedReadCounts),
    "\n"
)


############################################################
# BOOTSTRAP
############################################################

cat("Running bootstrap analysis...\n")

bootList <- BootList(
    geneNames,
    samplesNormalizedReadCounts,
    referenceNormalizedReadCounts,
    replicates = 10000
)

print(summary(bootList))


############################################################
# BACKGROUND NOISE
############################################################

cat("Estimating background noise...\n")

backgroundNoise <- Background(
    geneNames,
    samplesNormalizedReadCounts,
    referenceNormalizedReadCounts,
    bootList,
    replicates = 10000,
    significanceLevel = 0.001
)

print(summary(backgroundNoise))


############################################################
# REPORT
############################################################

cat("Generating CNV report...\n")

reportTables <- ReportTables(
    geneNames,
    samplesNormalizedReadCounts,
    referenceNormalizedReadCounts,
    bootList,
    backgroundNoise
)

print(reportTables)

############################################################
# CONVERT TO DATA FRAME
############################################################

cnvs <- as.data.frame(reportTables)

cat(
    "CNV table dimensions:",
    dim(cnvs),
    "\n"
)

print(names(cnvs))


############################################################
# FILTER PASSED CNVs
############################################################

if (!"Tumor.Passed" %in% names(cnvs)) {
    stop(
        "ERROR: Tumor.Passed column was not found in CNVPanelizer output."
    )
}

cnvs <- cnvs[cnvs$Tumor.Passed != 0, ]

cnvs$Tumor.Gene <- rownames(cnvs)


############################################################
# OUTPUT COLUMNS
############################################################

if (ncol(cnvs) >= 14) {

    cnvs <- cnvs[, c(14, 1, 11, 12)]

} else {

    warning(
        "Expected at least 14 columns in CNVPanelizer output. ",
        "Returning all available columns."
    )

}


############################################################
# OUTPUT
############################################################

file <- samplePrefix

fname <- gsub("_T$", "", file)

outputFile <- paste0(
    fname,
    ".copyNumberVariations.csv"
)

write.csv(
    cnvs,
    file = outputFile,
    row.names = FALSE
)


############################################################
# FINISHED
############################################################

cat("\n")
cat("====================================================\n")
cat("CNV ANALYSIS COMPLETED\n")
cat("Output:", outputFile, "\n")
cat("====================================================\n")
