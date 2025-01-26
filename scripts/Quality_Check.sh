#!/bin/bash

################################################################################
# Script: Quality_Check.sh
# Description: This script provides with quality control analysis on trimmed FASTQ 
#              files using FastQC and summarizes the results with MultiQC.
# Outputs: 
#   - Individual FastQC reports in the specified output directory.
#   - A summarized MultiQC report in the specified output directory.
# Author: Isra Boughanmi
# Date: 26/01/2025
################################################################################


# Making the script executable
chmod +x "$0"

# Defining paths
FASTQC_EXEC="/Applications/FastQC.app/Contents/MacOS/fastqc"  # Path to FastQC executable
TRIMMED_FASTQ_DIR="/Users/u5677580/Desktop/processed_fastq"   # Directory containing trimmed FASTQ files
FASTQC_OUTPUT_DIR="/Users/u5677580/Desktop/fastqc_results_trimmed"  # Directory for FastQC results
MULTIQC_OUTPUT_DIR="/Users/u5677580/Desktop/multiqc_trimmed_results"  # Directory for MultiQC results

# Creating necessary directories
mkdir -p "$FASTQC_OUTPUT_DIR"
mkdir -p "$MULTIQC_OUTPUT_DIR"

# Running FastQC on trimmed FASTQ files
echo "Running FastQC on trimmed FASTQ files..."
for FILE in "$TRIMMED_FASTQ_DIR"/*.fastq.gz; do
    echo "Processing $FILE with FastQC..."
    "$FASTQC_EXEC" "$FILE" --outdir="$FASTQC_OUTPUT_DIR"
done
echo "FastQC analysis completed. Reports saved in $FASTQC_OUTPUT_DIR."

# Running MultiQC on FastQC results
echo "Running MultiQC on FastQC results..."
multiqc "$FASTQC_OUTPUT_DIR" -o "$MULTIQC_OUTPUT_DIR"
if [ $? -eq 0 ]; then
    echo "MultiQC analysis completed. Reports saved in $MULTIQC_OUTPUT_DIR."
else
    echo "MultiQC analysis failed!"
    exit 1
fi
