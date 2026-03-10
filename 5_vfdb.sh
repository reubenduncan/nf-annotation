for genome in output_ysf_1/final_bins_copy/*.fa; do
    base=$(basename "$genome" .fa)

    abricate \
        --db vfdb \
        --minid 80 \
        --mincov 80 \
        "$genome" \
        > abricate_results/"${base}"_vfdb.tsv
done
