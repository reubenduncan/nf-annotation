#!/usr/bin/env bash
set -euo pipefail

INPUT_DIR=~/output_ysf_1/final_bins_copy
OUTPUT_DIR=~/6_mges

mkdir -p ~/6_mges

for genome in $INPUT_DIR/*.fa; do
    base=$(basename "$genome" .fa)

#    podman run \
#      --network host \
#      -v $INPUT_DIR:/data/input \
#      -v $OUTPUT_DIR:/data/output \ 
#      staphb/plasmidfinder \
#        -i /data/input/${base}.fa \
#        -o /data/output/${base} \
#        -p db/plasmidfinder_db

    mkdir $OUTPUT_DIR/$base

    plasmidfinder.py -i $genome -o $OUTPUT_DIR/$base -p db/plasmidfinder_db

done
