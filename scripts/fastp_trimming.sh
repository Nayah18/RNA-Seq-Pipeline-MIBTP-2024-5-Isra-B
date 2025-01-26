#!/bin/bash

################################################################################
# Script: fastp_trimming_pipeline.sh
# Description: This script performs adapter trimming and quality filtering on 
#              paired-end FASTQ files using fastp. It generates trimmed FASTQ 
#              files and detailed reports for each sample.
# Outputs: 
#   - Trimmed FASTQ files in the specified output directory.
#   - HTML and JSON reports for each sample in the report directory.
# Author: Isra Boughanmi
# Date: 26/01/2025
################################################################################


# Defining paths
FASTQ_DIR="/Users/u5677580/Downloads/data_014_extracted" 
OUTPUT_DIR="/Users/u5677580/Desktop/processed_fastq"     
REPORT_DIR="/Users/u5677580/Desktop/fastp_reports"       
THREADS=4                                                

# Creating output and report directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$REPORT_DIR"

# Looping through paired-end FASTQ files
for SAMPLE in $(ls $FASTQ_DIR/*_1.fastq.gz | sed 's/_1.fastq.gz//'); do
    BASENAME=$(basename "$SAMPLE")
    
    echo "Processing sample: $BASENAME"
    
    # Running fastp for paired-end trimming
    fastp \
        -i "${SAMPLE}_1.fastq.gz" \
        -I "${SAMPLE}_2.fastq.gz" \
        -o "$OUTPUT_DIR/${BASENAME}_1_trimmed.fastq.gz" \
        -O "$OUTPUT_DIR/${BASENAME}_2_trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread $THREADS \
        --html "$REPORT_DIR/${BASENAME}_fastp.html" \
        --json "$REPORT_DIR/${BASENAME}_fastp.json"
    
    echo "Finished processing $BASENAME" # To ensure all files have been processed
done

