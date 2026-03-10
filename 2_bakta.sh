for genome in output_ysf_1/final_bins/*.fa.gz; do
    base=$(basename "$genome" .fa.gz)

    bakta \
        --db ~/db/db \
        --output bakta_results/"$base" \
        --threads 8 \
        "$genome"
done
