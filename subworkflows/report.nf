// ─────────────────────────────────────────────────────────────────────────────
// subworkflows/report.nf
//
// GENERATE_REPORT  –  Read combined_results.parquet and produce a self-contained
// HTML report in the same dark-theme style as demo_report.html.
// No LLM interpretation — data, charts, and tables only.
//
// The Python logic lives in bin/generate_report.py (Nextflow adds bin/ to PATH
// automatically, so the script is called directly by name).
// ─────────────────────────────────────────────────────────────────────────────

process GENERATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    container  'quay.io/biocontainers/pandas:1.5.2'
    conda      'conda-forge::python=3.11 conda-forge::pandas conda-forge::pyarrow'

    input:
    path parquet

    output:
    path 'annotation_report.html', emit: html

    script:
    """
    generate_report.py ${parquet}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// REPORT workflow
// ─────────────────────────────────────────────────────────────────────────────
workflow REPORT {

    take:
    parquet

    main:
    GENERATE_REPORT(parquet)

    emit:
    html = GENERATE_REPORT.out.html
}
