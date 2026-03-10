#!/usr/bin/env bash
set -euo pipefail

mkdir -p amr_results

for genome in output_ysf_1/final_bins_copy/*.fa; do
    base=$(basename "$genome" .fa)

    

    amrfinder \
        -n "$genome" \
        -o amr_results/"$base"_amr.tsv \
        --threads 8
done

