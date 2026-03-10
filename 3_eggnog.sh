#!/usr/bin/env bash
set -euo pipefail

for proteome in bakta_results/*/*.faa; do
    base=$(basename "$proteome" .faa)

    emapper.py \
        -i "$proteome" \
        -o "$base" \
        --output_dir eggnog_results \
        --cpu 8
done
