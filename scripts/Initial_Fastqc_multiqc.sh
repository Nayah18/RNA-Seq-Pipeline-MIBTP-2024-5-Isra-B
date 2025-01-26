#!/bin/bash

################################################################################
# Script: Initial_Fastqc_Multiqc.sh
# Description: This script performs quality control on raw FASTQ files using 
#              FastQC and aggregates the results with MultiQC. It also extracts 
#              input data from a ZIP file if necessary for reproducibility.
# Outputs: 
#   - FastQC reports for individual FASTQ files.
#   - MultiQC report aggregating FastQC results.
# Author: Isra Boughanmi
# Date: 26/01/2025
################################################################################

# Making the script executable (if not already)
chmod +x "$0"

# Defining paths
FASTQC_EXEC="/Applications/FastQC.app/Contents/MacOS/fastqc"  # Path to FastQC
ZIP_FILE="/Users/u5677580/Downloads/data_014.zip"  # Path to the main zip file
EXTRACT_DIR="/Users/u5677580/Downloads/data_014_extracted"  # Directory to extract files
OUTPUT_DIR="/Users/u5677580/Desktop/fastqc_results"  # Directory to save FastQC results
MULTIQC_OUTPUT_DIR="/Users/u5677580/Desktop/multiqc_results"  # Directory for MultiQC report

# Creating necessary directories
mkdir -p "$EXTRACT_DIR"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$MULTIQC_OUTPUT_DIR"

# Extracting the main zip file with overwrite enabled
echo "Extracting $ZIP_FILE..."
if unzip -o "$ZIP_FILE" -d "$EXTRACT_DIR"; then
    echo "Extraction successful!"
else
    echo "Failed to extract $ZIP_FILE"
    exit 1
fi

# Navigating to the extracted directory
cd "$EXTRACT_DIR" || { echo "Extracted directory not found!"; exit 1; }

# Running FastQC for each FASTQ file
echo "Running FastQC for each FASTQ file..."
for FILE in *.fastq.gz; do
    echo "Processing $FILE with FastQC..."
    "$FASTQC_EXEC" "$FILE" --outdir="$OUTPUT_DIR"
done

echo "FastQC analysis completed. Reports saved in $OUTPUT_DIR."

# Running MultiQC to aggregate FastQC results
echo "Running MultiQC on FastQC results..."
multiqc "$OUTPUT_DIR" -o "$MULTIQC_OUTPUT_DIR"
if [ $? -eq 0 ]; then
    echo "MultiQC analysis completed. Aggregated report saved in $MULTIQC_OUTPUT_DIR."
else
    echo "MultiQC analysis failed!"
    exit 1
fi

echo "All quality control steps completed successfully!"
