#!/usr/bin/env bash
set -euo pipefail

KOFAM_DB=~/db/kofam
INPUT_DIR=~/output_ysf_1/NFAnnotation/1_bakta
OUTPUT_DIR=~/5_kofam

mkdir -p ~/5_kofam

for proteome in ~/output_ysf_1/NFAnnotation/1_bakta/*/*.fa.faa; do
    base=$(basename "$proteome" .fa.faa)

    echo "Annotating $base"

    podman run \
      --network host \
      -v $KOFAM_DB:/data/kofam_db \
      -v $INPUT_DIR:/data/input \
      -v $OUTPUT_DIR:/data/output \
      reubenduncan/kofam_scan:amd64 \
        -o /data/output/"$base"_kofam.tsv \
        -p /data/kofam_db/profiles \
        -k /data/kofam_db/ko_list \
        --cpu 8 \
        "/data/input/${base}/${base}.fa.faa"
done
