#!/bin/bash
# ==============================================
# Script: align_sort.sh
# Author: Isra Boughanmi
# Date: 26/01/2025
# Description:
#   This script aligns paired-end FASTQ files to a 
#   reference genome using HISAT2, converts SAM files 
#   to sorted BAM files, and generates alignment summaries.
# Requirements:
#   - HISAT2
#   - SAMtools
#   - Prebuilt HISAT2 reference index
# ==============================================

# Making the script executable (if not already)
chmod +x "$0"

# Paths
REFERENCE_DIR="/Users/u5677580/RNAseq_human/reference"
INDEX_PREFIX="/Users/u5677580/RNAseq_human/reference/grch37_tran/genome_tran"
FASTQ_DIR="/Users/u5677580/Downloads/data_014_extracted"
ALIGNED_DIR="/Users/u5677580/Desktop/aligned"
SORTED_BAM_DIR="/Users/u5677580/Desktop/sorted_bam"
SUMMARY_DIR="/Users/u5677580/Desktop/alignment_summaries"

# Creating necessary directories
echo "Creating directories for aligned SAM, sorted BAM files, and alignment summaries..."
mkdir -p "$ALIGNED_DIR" "$SORTED_BAM_DIR" "$SUMMARY_DIR"

# Aligning raw FASTQ files using HISAT2
echo "Aligning raw FASTQ files using HISAT2..."
for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    BASENAME=$(basename "$R1" _1.fastq.gz)
    R2="$FASTQ_DIR/${BASENAME}_2.fastq.gz"

    # Checking if paired-end files exist (was necessary for debugging initially, but not required to run the script)
    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "Processing sample: $BASENAME"

        hisat2 -x "$INDEX_PREFIX" \
               -1 "$R1" \
               -2 "$R2" \
               -S "$ALIGNED_DIR/${BASENAME}.sam" \
               --threads 4 \
               --summary-file "$SUMMARY_DIR/${BASENAME}_alignment_summary.txt"

        echo "Finished aligning $BASENAME"
    else
        echo "Warning: Paired files for $BASENAME not found. Skipping..."
    fi
done

# Converting and sorting SAM files to BAM format
echo "Converting and sorting SAM files to BAM format..."
for SAM_FILE in "$ALIGNED_DIR"/*.sam; do
    BASENAME=$(basename "$SAM_FILE" .sam)
    echo "Sorting and converting $SAM_FILE"

    samtools sort -@ 4 -o "$SORTED_BAM_DIR/${BASENAME}_sorted.bam" "$SAM_FILE"
    samtools index "$SORTED_BAM_DIR/${BASENAME}_sorted.bam"

    echo "Finished processing $BASENAME"
done
