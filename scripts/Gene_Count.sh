# ==============================================
# Script: Gene_Count.sh
# Author: Isra Boughanmi
# Date: 26/01/2025
# Description:
#   This script uses LiBiNorm to calculate gene counts 
#   from sorted BAM files, based on a specified GTF file. 
#   The results are saved in a specified directory.
# Requirements:
#   - LiBiNorm
#   - Sorted BAM files
#   - GTF file
# ==============================================
#!/bin/bash

# Making the script executable (if not already)
chmod +x "$0"

# Paths
SORTED_BAM_DIR="/Users/u5677580/Desktop/sorted_bam"       # Directory containing sorted BAM files
REFERENCE_DIR="/Users/u5677580/RNAseq_human/reference"    # Path to reference directory
GTF_FILE="$REFERENCE_DIR/homo_sapiens.grch37.75.gtf"      # Path to the GTF file
COUNTS_DIR="/Users/u5677580/Desktop/gene_counts"          # Directory to save gene count files

# Creating output directory for counts
echo "Creating directory for gene counts..."
mkdir -p "$COUNTS_DIR"

# Running LiBiNorm to generate gene counts, with debugging messages
echo "Generating gene counts using LiBiNorm..."
for BAM_FILE in "$SORTED_BAM_DIR"/*.bam; do
    if [ -e "$BAM_FILE" ]; then
        # Remove the "_sorted.bam" suffix to match the count file naming convention
        BASENAME=$(basename "$BAM_FILE" _sorted.bam)
        OUTPUT_FILE="$COUNTS_DIR/${BASENAME}_counts.txt"

        echo "Processing $BAM_FILE for gene quantification..."

        LiBiNorm count -z -r pos -i gene_name -s reverse \
            "$BAM_FILE" "$GTF_FILE" > "$OUTPUT_FILE"

        echo "Gene counts saved to $OUTPUT_FILE"
    else
        echo "No BAM files found in $SORTED_BAM_DIR. Skipping..."
        exit 1
    fi
done

