// ─────────────────────────────────────────────────────────────────────────────
// subworkflows/merge.nf
//
// MERGE_RESULTS  –  Combine per-sample TSVs from all annotation tools into a
// single Parquet file.  Uses filename glob patterns so any number of samples
// can be collected into each input channel without the space-separated
// filename expansion problem.
// ─────────────────────────────────────────────────────────────────────────────

process MERGE_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'
    container 'python:3.11-slim'
    conda     'conda-forge::python=3.11 conda-forge::pandas conda-forge::pyarrow'
    when: params.run_merge

    input:
    path amrfinder_tsvs
    path kofam_tsvs
    path card_rgi_tsvs
    path plasmidfinder_tsvs
    path integronfinder_tsvs
    path mob_suite_tsvs
    path isescan_tsvs

    output:
    path 'combined_results.parquet', emit: parquet

    script:
    """
    pip install -q pandas pyarrow

    python3 << 'PYEOF'
import pandas as pd
import glob
import sys

# Map filename pattern -> source_tool label.
# All per-sample TSVs are staged into the working directory by Nextflow,
# so glob patterns reliably match across any number of samples.
TOOL_PATTERNS = [
    ('*_amr.tsv',          'amrfinder'),
    ('*_kofam.tsv',        'kofam'),
    ('*_rgi.tsv',          'card_rgi'),
    ('*_plasmidfinder.tsv','plasmidfinder'),
    ('*_integrons.tsv',    'integronfinder'),
    ('*_mob.tsv',          'mob_suite'),
    ('*_isescan.tsv',      'isescan'),
]

def read_tsv(path):
    return pd.read_csv(path, sep='\\t', comment='#', quotechar='"', quoting=1)

all_frames = []
for pattern, tool in TOOL_PATTERNS:
    for f in glob.glob(pattern):
        try:
            d = read_tsv(f)
            d['source_tool'] = tool
            all_frames.append(d)
            print(f"  {tool}: {f} ({len(d)} rows)", file=sys.stderr)
        except Exception as e:
            print(f"Warning – could not read {f}: {e}", file=sys.stderr)

if all_frames:
    combined = pd.concat(all_frames, axis=0, ignore_index=True, sort=True)
    cols = combined.columns.tolist()
    cols.remove('source_tool')
    combined = combined[['source_tool'] + cols]
    combined.to_parquet('combined_results.parquet', index=False)
    print(f"Merged {len(all_frames)} file(s) → {len(combined):,} rows total", file=sys.stderr)
else:
    print("Warning: no TSV files found to merge", file=sys.stderr)
    pd.DataFrame({'source_tool': pd.Series([], dtype='str')}).to_parquet(
        'combined_results.parquet', index=False)
PYEOF
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// MERGE workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow MERGE {

    take:
    amrfinder_tsvs
    kofam_tsvs
    card_rgi_tsvs
    plasmidfinder_tsvs
    integronfinder_tsvs
    mob_suite_tsvs
    isescan_tsvs

    main:
    MERGE_RESULTS(
        amrfinder_tsvs,
        kofam_tsvs,
        card_rgi_tsvs,
        plasmidfinder_tsvs,
        integronfinder_tsvs,
        mob_suite_tsvs,
        isescan_tsvs
    )

    emit:
    parquet = MERGE_RESULTS.out.parquet
}
