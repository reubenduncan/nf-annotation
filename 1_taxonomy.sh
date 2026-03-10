#!/bin/bash

# Script to filter bins by quality and run gtdbtk classify_wf
set -e

QUALITY_REPORT="/gpfs/home/s441125/output/quality_report.tsv"
BINS_DIR="/gpfs/home/s441125/output/final_bins"
WORK_DIR="/gpfs/home/s441125/output/gtdbtk_work"
OUTPUT_DIR="/gpfs/home/s441125/output/gtdbtk_out"

echo "Creating work directory..."
mkdir -p "$WORK_DIR"

echo "Filtering genomes with >= 80% completeness and <= 5% contamination..."
# Skip header line and filter based on criteria
awk -F'\t' 'NR > 1 && $2 >= 80 && $3 <= 5 {print $1}' "$QUALITY_REPORT" | while read genome_name; do
    # Find the corresponding gzipped file
    gz_file="$BINS_DIR/${genome_name}.fa.gz"
    
    if [ -f "$gz_file" ]; then
        echo "Processing: $genome_name"
        # Unzip directly to work directory
        gunzip -c "$gz_file" > "$WORK_DIR/${genome_name}.fa"
    else
        echo "Warning: File not found for $genome_name"
    fi
done

echo "Running gtdbtk classify_wf..."
gtdbtk classify_wf \
  --genome_dir "$WORK_DIR" \
  --out_dir "$OUTPUT_DIR" \
  --cpus 16 \
  --extension fa

echo "Cleaning up work directory..."
rm -rf "$WORK_DIR"

echo "Done! GTDB-Tk results saved to: $OUTPUT_DIR"
