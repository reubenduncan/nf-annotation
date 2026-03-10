#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================
// Parameters
// ============================================================
params.input              = null          // Input dir: flat *.fa.gz files, or subdirs each containing *.fa.gz files
params.quality_report     = null          // (optional) TSV: col1=name, col2=completeness, col3=contamination
params.bakta_db           = null          // Path to Bakta database directory
params.gtdbtk_db          = null          // Path to GTDB-Tk reference data directory
params.kofam_db           = null          // Path to KofamScan db dir (profiles/ + ko_list)
params.eggnog_db          = null          // Path to EggNOG-mapper data directory
params.plasmidfinder_db   = '/plasmidfinder_db' // PlasmidFinder db (default: in container)
params.db_dir             = 'databases'         // Where auto-downloaded databases are stored
params.outdir             = 'results'
params.min_completeness   = 80
params.max_contamination  = 5

// Pipeline stages — set to true to run each stage
params.run_taxonomy   = true
params.run_bakta      = true
params.run_eggnog     = true
params.run_amrfinder  = true
params.run_vfdb       = true
params.run_kofam      = true
params.run_mge        = true
params.run_antismash  = true
params.run_merge      = true   // final merge to Parquet

// ============================================================
// Validate required parameters
// Databases (bakta_db, gtdbtk_db, kofam_db, eggnog_db) are optional:
// if not provided they will be auto-downloaded to params.db_dir.
// ============================================================
def required_params = ['input']
required_params.each { p ->
    if (!params[p]) error "Missing required parameter: --${p}"
}

// ============================================================
// DECOMPRESS
// Gunzip each .fa.gz bin into a plain .fa file for downstream use.
// ============================================================
process DECOMPRESS {
    tag "$name"
    container 'ubuntu:24.04'

    input:
    tuple val(name), path(gz_file)

    output:
    tuple val(name), path("${name}.fa")

    script:
    """
    gunzip -c ${gz_file} > ${name}.fa
    """
}

// ============================================================
// CONCAT_BIN
// Subfolder input mode: concatenate all *.fa.gz files within a subdirectory
// into a single .fa.gz, which is then passed to DECOMPRESS like any other bin.
// Concatenated gzip streams are valid: gunzip handles multiple streams correctly.
// ============================================================
process CONCAT_BIN {
    tag "$name"
    container 'ubuntu:24.04'

    input:
    tuple val(name), path(gz_files)

    output:
    tuple val(name), path("${name}.fa.gz")

    script:
    """
    cat ${gz_files} > ${name}.fa.gz
    """
}

// ============================================================
// DATABASE DOWNLOADS
// Each process is only invoked when the user has not supplied
// a pre-built database path via the corresponding parameter.
// storeDir caches the result so subsequent runs skip the download.
// ============================================================

process DOWNLOAD_BAKTA_DB {
    storeDir "${params.db_dir}/bakta"
    container 'oschwengers/bakta:latest'
    label 'medium_cpu'

    output:
    path 'db', emit: db

    script:
    """
    bakta_db download --output db --type full
    """
}

process DOWNLOAD_GTDBTK_DB {
    storeDir "${params.db_dir}/gtdbtk"
    container 'ecogenomics/gtdbtk:2.4.0'
    label 'high_cpu'

    output:
    path 'gtdbtk_db', emit: db

    script:
    // download-db.sh is bundled in the ecogenomics/gtdbtk image;
    // it accepts a target directory as its only argument.
    // Warning: the GTDB-Tk reference package is ~66 GB.
    """
    download-db.sh gtdbtk_db
    """
}

process DOWNLOAD_KOFAM_DB {
    storeDir "${params.db_dir}/kofam"
    container 'ubuntu:24.04'

    output:
    path 'kofam_db', emit: db

    script:
    """
    mkdir -p kofam_db
    wget -q ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
    wget -q ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
    tar -xzf profiles.tar.gz
    mv profiles kofam_db/
    gunzip ko_list.gz
    mv ko_list kofam_db/
    rm profiles.tar.gz
    """
}

process DOWNLOAD_EGGNOG_DB {
    storeDir "${params.db_dir}/eggnog"
    container 'nanozoo/eggnog-mapper:2.1.13--c16a7d2'

    output:
    path 'eggnog_db', emit: db

    script:
    """
    mkdir -p eggnog_db
    download_eggnog_data.py --data_dir eggnog_db -y
    """
}

// ============================================================
// TAXONOMY  (Step 1)
// Collect all quality-filtered genomes and classify with GTDB-Tk.
// Runs as a single job over all genomes (GTDB-Tk requires a directory).
// Container: ecogenomics/gtdbtk  https://hub.docker.com/r/ecogenomics/gtdbtk
// ============================================================
process TAXONOMY {
    tag "gtdbtk"
    label 'high_cpu'
    publishDir "${params.outdir}/taxonomy", mode: 'copy'
    container 'nanozoo/gtdbtk:2.4.0--02c00d5'
    when: params.run_taxonomy

    input:
    path fa_files       // collected list of .fa files
    path gtdbtk_db

    output:
    path 'gtdbtk_out/'

    script:
    """
    export GTDBTK_DATA_PATH=${gtdbtk_db}
    mkdir -p genome_input
    cp ${fa_files} genome_input/

    gtdbtk classify_wf \\
        --genome_dir genome_input \\
        --out_dir    gtdbtk_out  \\
        --cpus       ${task.cpus} \\
        --extension  fa
    """
}

// ============================================================
// BAKTA  (Step 2)
// Structural and functional annotation of each genome.
// Container: oschwengers/bakta  https://hub.docker.com/r/oschwengers/bakta
// ============================================================
process BAKTA {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/bakta", mode: 'copy'
    container 'oschwengers/bakta:latest'
    when: params.run_bakta

    input:
    tuple val(name), path(fa_file)
    path bakta_db

    output:
    tuple val(name), path("${name}/"),          emit: dir
    tuple val(name), path("${name}/${name}.faa"), emit: faa

    script:
    """
    bakta \\
        --db      ${bakta_db} \\
        --output  ${name}     \\
        --prefix  ${name}     \\
        --threads ${task.cpus} \\
        ${fa_file}
    """
}

// ============================================================
// EGGNOG  (Step 3)
// Orthology and functional annotation via EggNOG-mapper.
// Depends on Bakta .faa protein output.
// Container: nanozoo/eggnog-mapper  https://hub.docker.com/r/nanozoo/eggnog-mapper
// ============================================================
process EGGNOG {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/eggnog", mode: 'copy'
    container 'nanozoo/eggnog-mapper:2.1.13--c16a7d2'
    when: params.run_eggnog && params.run_bakta

    input:
    tuple val(name), path(faa_file)
    path eggnog_db

    output:
    tuple val(name), path("${name}.*")

    script:
    """
    emapper.py \\
        -i          ${faa_file} \\
        -o          ${name}     \\
        --output_dir . \\
        --data_dir  ${eggnog_db} \\
        --cpu       ${task.cpus}
    """
}

// ============================================================
// AMRFINDER  (Step 4)
// AMR gene detection with NCBI AMRFinderPlus.
// Container: staphb/ncbi-amrfinderplus  https://hub.docker.com/r/staphb/ncbi-amrfinderplus
// ============================================================
process AMRFINDER {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/amrfinder", mode: 'copy'
    container 'staphb/ncbi-amrfinderplus:latest'
    when: params.run_amrfinder

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}_amr.tsv")

    script:
    """
    amrfinder \\
        -n ${fa_file}            \\
        -o ${name}_amr.tsv \\
        --threads ${task.cpus}
    """
}

// ============================================================
// VFDB  (Step 5)
// Virulence factor annotation via abricate against VFDB.
// Container: staphb/abricate  https://hub.docker.com/r/staphb/abricate
// ============================================================
process VFDB {
    tag "$name"
    publishDir "${params.outdir}/vfdb", mode: 'copy'
    container 'staphb/abricate:latest'
    when: params.run_vfdb

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}_vfdb.tsv")

    script:
    """
    abricate \\
        --db     vfdb \\
        --minid  80   \\
        --mincov 80   \\
        ${fa_file}    \\
        > ${name}_vfdb.tsv
    """
}

// ============================================================
// KOFAM  (Step 6)
// KEGG Orthology annotation via KofamScan.
// Outputs space-separated format; KOFAM_REFORMAT will convert to proper TSV.
// Depends on Bakta .faa protein output.
// Container: reubenduncan/kofam_scan  https://hub.docker.com/r/reubenduncan/kofam_scan
// ============================================================
process KOFAM {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/kofam_raw", mode: 'copy'
    container 'reubenduncan/kofam_scan:amd64'
    when: params.run_kofam && params.run_bakta

    input:
    tuple val(name), path(faa_file)
    path kofam_db

    output:
    tuple val(name), path("${name}_kofam.tsv"), emit: raw

    script:
    """
    exec_annotation \\
        -o ${name}_kofam.tsv     \\
        -p ${kofam_db}/profiles  \\
        -k ${kofam_db}/ko_list   \\
        --cpu ${task.cpus}       \\
        ${faa_file}
    """
}

// ============================================================
// KOFAM_REFORMAT
// Convert KofamScan space-separated output to proper TSV.
// Input has format: gene_name  KO  threshold  score  E-value  KO_definition
// Output: tab-separated with proper header.
// ============================================================
process KOFAM_REFORMAT {
    tag "$name"
    publishDir "${params.outdir}/kofam", mode: 'copy'
    container 'python:3.11-slim'
    when: params.run_kofam

    input:
    tuple val(name), path(kofam_raw)

    output:
    tuple val(name), path("${name}_kofam.tsv")

    script:
    """
    python3 << 'EOF'
import sys

# Read the space-separated file, skip comment lines
lines = []
with open('${kofam_raw}', 'r') as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        lines.append(line)

# Write proper TSV: split on whitespace (multiple spaces/tabs = single delimiter)
with open('${name}_kofam.tsv', 'w') as out:
    # Header
    out.write('gene_name\tKO\tthreshold\tscore\tE_value\tKO_definition\\n')

    for line in lines:
        fields = line.split()
        if len(fields) >= 6:
            # Take first 5 columns as-is, join remaining as KO definition in case it has spaces
            gene_name = fields[0]
            ko = fields[1]
            threshold = fields[2]
            score = fields[3]
            e_value = fields[4]
            ko_def = ' '.join(fields[5:])

            out.write(f"{gene_name}\\t{ko}\\t{threshold}\\t{score}\\t{e_value}\\t{ko_def}\\n")
EOF
    """
}

// ============================================================
// MGE  (Step 7)
// Mobile genetic element / plasmid detection via PlasmidFinder.
// Outputs JSON results; MGE_TO_TSV will convert to TSV.
// Container: staphb/plasmidfinder  https://hub.docker.com/r/staphb/plasmidfinder
// Database bundled in container at /plasmidfinder_db by default;
// override with --plasmidfinder_db to use an external copy.
// ============================================================
process MGE {
    tag "$name"
    publishDir "${params.outdir}/mge_raw", mode: 'copy'
    container 'staphb/plasmidfinder:latest'
    when: params.run_mge

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}/results.json"), emit: json

    script:
    """
    mkdir -p ${name}
    plasmidfinder.py \\
        -i ${fa_file}                      \\
        -o ${name}                         \\
        -p ${params.plasmidfinder_db}
    """
}

// ============================================================
// MGE_TO_TSV
// Convert PlasmidFinder JSON output to TSV format.
// Schema: genome, contig, match_id, match_name, coverage, identity, hit_length
// ============================================================
process MGE_TO_TSV {
    tag "$name"
    publishDir "${params.outdir}/mge", mode: 'copy'
    container 'python:3.11-slim'
    when: params.run_mge

    input:
    tuple val(name), path(results_json)

    output:
    tuple val(name), path("${name}_mge.tsv")

    script:
    """
    python3 << 'EOF'
import json
import sys

# Load PlasmidFinder JSON
with open('${results_json}', 'r') as f:
    data = json.load(f)

# Write TSV with consistent schema
with open('${name}_mge.tsv', 'w') as out:
    out.write('genome\\tcontig\\tmatch_id\\tmatch_name\\tcoverage\\tidentity\\thit_length\\n')

    # PlasmidFinder JSON structure: "results" -> contig_name -> list of matches
    results = data.get('results', {})
    for contig_name, matches in results.items():
        if isinstance(matches, list):
            for match in matches:
                genome = '${name}'
                contig = contig_name
                match_id = match.get('match_id', '')
                match_name = match.get('match_name', '')
                coverage = match.get('coverage', '')
                identity = match.get('identity', '')
                hit_length = match.get('hit_length', '')

                out.write(f"{genome}\\t{contig}\\t{match_id}\\t{match_name}\\t{coverage}\\t{identity}\\t{hit_length}\\n")
EOF
    """
}

// ============================================================
// ANTISMASH  (Step 8)
// Secondary metabolite biosynthetic gene cluster prediction.
// Container: antismash/antismash  https://hub.docker.com/r/antismash/antismash
// ============================================================
process ANTISMASH {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/antismash", mode: 'copy'
    container 'antismash/antismash:latest'
    when: params.run_antismash

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}/")

    script:
    """
    antismash \\
        ${fa_file}                     \\
        --output-dir     ${name}       \\
        --genefinding-tool prodigal    \\
        --cpus           ${task.cpus}
    """
}

// ============================================================
// MERGE_RESULTS
// Merge all TSV outputs into a single Parquet file.
// Handles ## comment lines; adds source_tool column.
// ============================================================
process MERGE_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'
    container 'python:3.11-slim'
    when: params.run_merge

    input:
    path amrfinder_tsvs
    path vfdb_tsvs
    path kofam_tsvs
    path mge_tsvs

    output:
    path 'combined_results.parquet'

    script:
    """
    python3 << 'EOF'
import pandas as pd
import glob
import os

# Collect all TSVs with their source tool labels
all_data = []

# Helper to read TSVs safely (skip ## comment lines)
def read_tsv_safe(file_path):
    df = pd.read_csv(file_path, sep='\\t', comment='#', quotechar='"', quoting=1)
    return df

# Read AMRFinder results
for f in glob.glob('${amrfinder_tsvs}'):
    if f.endswith('.tsv'):
        try:
            df = read_tsv_safe(f)
            df['source_tool'] = 'amrfinder'
            all_data.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=__import__('sys').stderr)

# Read VFDB results
for f in glob.glob('${vfdb_tsvs}'):
    if f.endswith('.tsv'):
        try:
            df = read_tsv_safe(f)
            df['source_tool'] = 'vfdb'
            all_data.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=__import__('sys').stderr)

# Read KOFAM results
for f in glob.glob('${kofam_tsvs}'):
    if f.endswith('.tsv'):
        try:
            df = read_tsv_safe(f)
            df['source_tool'] = 'kofam'
            all_data.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=__import__('sys').stderr)

# Read MGE (PlasmidFinder) results
for f in glob.glob('${mge_tsvs}'):
    if f.endswith('.tsv'):
        try:
            df = read_tsv_safe(f)
            df['source_tool'] = 'mge_plasmidfinder'
            all_data.append(df)
        except Exception as e:
            print(f"Warning: Could not read {f}: {e}", file=__import__('sys').stderr)

if all_data:
    # Concatenate with alignment='minimal' (aligns columns by name only)
    combined = pd.concat(all_data, axis=0, ignore_index=True, sort=True)
    # Move source_tool to the front
    cols = combined.columns.tolist()
    cols.remove('source_tool')
    combined = combined[['source_tool'] + cols]

    # Export to Parquet
    combined.to_parquet('combined_results.parquet', index=False)
    print(f"Merged {len(all_data)} TSV files into combined_results.parquet ({len(combined)} rows)")
else:
    print("Warning: No TSV files found to merge")
    # Create empty parquet
    pd.DataFrame({'source_tool': []}).to_parquet('combined_results.parquet', index=False)
EOF
    """
}

// ============================================================
// Workflow
// ============================================================
workflow {

    // Build bins channel.
    // Flat mode:     input/ contains *.fa.gz directly → each file is one sample.
    // Subfolder mode: input/ contains subdirs with *.fa.gz → each subdir is one
    //                 sample; files within it are concatenated before processing.
    def input_path = file(params.input)
    def has_subdirs = input_path.listFiles()?.any { it.isDirectory() } ?: false

    def bins_ch
    if (has_subdirs) {
        def raw_ch = Channel
            .fromPath("${params.input}/*/*.fa.gz")
            .map { f -> [f.parent.name, f] }
            .groupTuple()
        CONCAT_BIN(raw_ch)
        bins_ch = CONCAT_BIN.out
    } else {
        bins_ch = Channel
            .fromPath("${params.input}/*.fa.gz")
            .map { f -> [f.name.replace('.fa.gz', ''), f] }
    }

    // Optionally filter bins by quality report; if not supplied all bins proceed.
    def filtered_bins
    if (params.quality_report) {
        def passing_names = Channel
            .fromPath(params.quality_report)
            .splitCsv(sep: '\t', skip: 1)
            .filter { row ->
                row[1].toFloat() >= params.min_completeness.toFloat() &&
                row[2].toFloat() <= params.max_contamination.toFloat()
            }
            .map { row -> [row[0], true] }
        filtered_bins = bins_ch
            .join(passing_names)
            .map { name, gz, _flag -> [name, gz] }
    } else {
        filtered_bins = bins_ch
    }

    // Decompress each filtered bin once; reuse fa_ch for all downstream steps
    DECOMPRESS(filtered_bins)
    fa_ch = DECOMPRESS.out   // (name, fa_file)

    // --- Database channels ---
    // Use user-supplied paths when provided, otherwise auto-download.
    def bakta_db_ch
    if (params.bakta_db) {
        bakta_db_ch = Channel.value(file(params.bakta_db))
    } else {
        DOWNLOAD_BAKTA_DB()
        bakta_db_ch = DOWNLOAD_BAKTA_DB.out.db.first()
    }

    def gtdbtk_db_ch
    if (params.gtdbtk_db) {
        gtdbtk_db_ch = Channel.value(file(params.gtdbtk_db))
    } else {
        DOWNLOAD_GTDBTK_DB()
        gtdbtk_db_ch = DOWNLOAD_GTDBTK_DB.out.db.first()
    }

    def kofam_db_ch
    if (params.kofam_db) {
        kofam_db_ch = Channel.value(file(params.kofam_db))
    } else {
        DOWNLOAD_KOFAM_DB()
        kofam_db_ch = DOWNLOAD_KOFAM_DB.out.db.first()
    }

    def eggnog_db_ch
    if (params.eggnog_db) {
        eggnog_db_ch = Channel.value(file(params.eggnog_db))
    } else {
        DOWNLOAD_EGGNOG_DB()
        eggnog_db_ch = DOWNLOAD_EGGNOG_DB.out.db.first()
    }

    // Step 1 — Taxonomy: collect all FAs into one GTDB-Tk classify_wf run
    TAXONOMY(
        fa_ch.map { name, fa -> fa }.collect(),
        gtdbtk_db_ch
    )

    // Step 2 — Bakta: per-genome structural + functional annotation
    BAKTA(fa_ch, bakta_db_ch)

    // Steps 3 & 6 — depend on Bakta protein sequences
    EGGNOG(BAKTA.out.faa, eggnog_db_ch)
    KOFAM(BAKTA.out.faa,  kofam_db_ch)

    // Step 6b — Reformat KOFAM output from space-separated to proper TSV
    KOFAM_REFORMAT(KOFAM.out.raw)

    // Steps 4, 5, 7, 8 — independent per-genome, run in parallel
    AMRFINDER(fa_ch)
    VFDB(fa_ch)
    MGE(fa_ch)
    ANTISMASH(fa_ch)

    // Step 7b — Convert MGE JSON to TSV
    MGE_TO_TSV(MGE.out.json)

    // Final step — merge all TSVs into a Parquet file
    if (params.run_merge) {
        // Collect all TSVs for merging
        amrfinder_tsvs = AMRFINDER.out.map { name, tsv -> tsv }.collect()
        vfdb_tsvs = VFDB.out.map { name, tsv -> tsv }.collect()
        kofam_tsvs = KOFAM_REFORMAT.out.map { name, tsv -> tsv }.collect()
        mge_tsvs = MGE_TO_TSV.out.map { name, tsv -> tsv }.collect()

        MERGE_RESULTS(amrfinder_tsvs, vfdb_tsvs, kofam_tsvs, mge_tsvs)
    }
}
