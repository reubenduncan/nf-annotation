#!/usr/bin/env bash
set -euo pipefail

mkdir -p 8_antismash

for genome in output_ysf_1/final_bins_copy/*.fa; do
    base=$(basename "$genome" .fa)

    antismash \
        "$genome" \
        --output-dir 8_antismash/"$base" \
        --genefinding-tool prodigal
done
