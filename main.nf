#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { COLLECT_INPUTS } from './subworkflows/input'
include { MERGE }         from './subworkflows/merge'
include { REPORT }        from './subworkflows/report'

// ============================================================
// Parameters
// ============================================================
params.input              = null          // Input dir — see subworkflows/input.nf for supported formats and layout modes
params.quality_report     = null          // (optional) TSV: col1=name, col2=completeness, col3=contamination
params.bakta_db           = null          // Path to Bakta database directory
params.gtdbtk_db          = null          // Path to GTDB-Tk reference data directory
params.kofam_db           = null          // Path to KofamScan db dir (profiles/ + ko_list)
params.eggnog_db          = null          // Path to EggNOG-mapper data directory
params.card_db            = null          // Path to CARD db dir (card.json + variants/); auto-downloaded if absent
params.plasmidfinder_db   = '/plasmidfinder_db' // PlasmidFinder db (default: in container)
params.db_dir             = 'databases'         // Where auto-downloaded databases are stored
params.outdir             = 'results'
params.min_completeness   = 80
params.max_contamination  = 5

// Pipeline stages — set to true to run each stage
params.run_taxonomy       = true
params.run_bakta          = true
params.run_eggnog         = true
params.run_amrfinder      = true
params.run_card_rgi       = true   // AMR via CARD RGI (replaces abricate/VFDB)
params.run_kofam          = true
params.run_metabolic      = true   // Metabolic pathway reconstruction (METABOLIC)
params.run_metabolishmm  = true   // HMM-based metabolic marker detection (metabolisHMM)
params.run_microtrait     = true   // Microbial functional trait inference (microTrait)
params.run_plasmidfinder  = true   // Plasmid detection (PlasmidFinder)
params.run_integronfinder = true   // Integron detection (IntegronFinder)
params.run_mob_suite      = true   // Plasmid reconstruction and MOB typing (mob_suite)
params.run_isescan        = true   // IS element detection (ISEScan)
params.run_antismash      = true
params.run_merge          = true   // final merge to Parquet
params.run_report         = true   // HTML annotation report

// ============================================================
// Validate required parameters
// Databases (bakta_db, gtdbtk_db, kofam_db, eggnog_db, card_db) are optional:
// if not provided they will be auto-downloaded to params.db_dir.
// ============================================================
def required_params = ['input']
required_params.each { p ->
    if (!params[p]) error "Missing required parameter: --${p}"
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
    conda     'bioconda::bakta'
    label 'medium_cpu'

    output:
    path 'db/db', emit: db

    script:
    """
    bakta_db download --output db --type full
    """
}

process DOWNLOAD_GTDBTK_DB {
    storeDir "${params.db_dir}/gtdbtk"
    container 'ecogenomics/gtdbtk:2.4.0'
    conda     'bioconda::gtdbtk'
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
    container 'reubenduncan/kofam_scan:amd64'
    conda     'bioconda::kofamscan conda-forge::wget'

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
    conda     'bioconda::eggnog-mapper'

    output:
    path 'eggnog_db', emit: db

    script:
    """
    mkdir -p eggnog_db
    download_eggnog_data.py --data_dir eggnog_db -y
    """
}

process DOWNLOAD_CARD_DB {
    storeDir "${params.db_dir}/card"
    container 'finlaymaguire/rgi:latest'
    conda     'bioconda::rgi conda-forge::wget'

    output:
    path 'card_db', emit: db

    script:
    """
    mkdir -p card_db
    cd card_db
    wget -q https://card.mcmaster.ca/latest/data -O card_data.tar.bz2
    tar -xjf card_data.tar.bz2
    rm card_data.tar.bz2
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
    conda     'bioconda::gtdbtk'
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
    conda     'bioconda::bakta'
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
    conda     'bioconda::eggnog-mapper'
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
    conda     'bioconda::ncbi-amrfinderplus'
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
// CARD_RGI  (Step 5)
// Comprehensive AMR annotation via CARD Resistance Gene Identifier.
// Depends on Bakta .faa protein output.
// Container: finlaymaguire/rgi  https://hub.docker.com/r/finlaymaguire/rgi
// ============================================================
process CARD_RGI {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/card_rgi", mode: 'copy'
    container 'finlaymaguire/rgi:latest'
    conda     'bioconda::rgi'
    when: params.run_card_rgi && params.run_bakta

    input:
    tuple val(name), path(faa_file)
    path card_db

    output:
    tuple val(name), path("${name}_rgi.tsv"), emit: tsv

    script:
    """
    rgi load --card_json ${card_db}/card.json --local
    rgi main \\
        --input_sequence ${faa_file} \\
        --output_file    ${name}_rgi \\
        --input_type     protein     \\
        --alignment_tool DIAMOND     \\
        --num_threads    ${task.cpus} \\
        --clean
    mv ${name}_rgi.txt ${name}_rgi.tsv
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
    conda     'bioconda::kofamscan'
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
    conda     'conda-forge::python=3.11'
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
// METABOLIC  (Step 7)
// Metabolic and biogeochemical functional annotation.
// Reconstructs metabolic pathways from genome nucleotide sequences.
// Container: jolespin/metabolic  https://hub.docker.com/r/jolespin/metabolic
// ============================================================
process METABOLIC {
    tag "$name"
    label 'high_cpu'
    publishDir "${params.outdir}/metabolic", mode: 'copy'
    container 'jolespin/metabolic:v4.0'
    conda     'bioconda::metabolic'
    when: params.run_metabolic

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}_metabolic/")

    script:
    """
    mkdir -p input_genomes
    cp ${fa_file} input_genomes/
    METABOLIC-G.pl \\
        -in-gn input_genomes \\
        -o     ${name}_metabolic \\
        -t     ${task.cpus}
    """
}

// ============================================================
// METABOLISHMM  (Step 8)
// HMM-based detection of metabolic pathway marker genes.
// Depends on Bakta .faa protein output.
// Container: elizabethmcd/metabolishmm  https://hub.docker.com/r/elizabethmcd/metabolishmm
// ============================================================
process METABOLISHMM {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/metabolishmm", mode: 'copy'
    container 'elizabethmcd/metabolishmm:latest'
    conda     'bioconda::metabolishmm'
    when: params.run_metabolishmm && params.run_bakta

    input:
    tuple val(name), path(faa_file)

    output:
    tuple val(name), path("${name}_metabolishmm/")

    script:
    """
    metabolishmm run-metabolic \\
        --input   ${faa_file}           \\
        --output  ${name}_metabolishmm \\
        --threads ${task.cpus}
    """
}

// ============================================================
// MICROTRAIT  (Step 9)
// Microbial functional trait inference from genome sequences.
// Container: ukaraoz/microtrait  https://hub.docker.com/r/ukaraoz/microtrait
// ============================================================
process MICROTRAIT {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/microtrait", mode: 'copy'
    container 'ukaraoz/microtrait:latest'
    conda     'bioconda::r-microtrait conda-forge::r-base=4'
    when: params.run_microtrait

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}_microtrait/")

    script:
    """
    microtrait infer \\
        --input  ${fa_file}          \\
        --outdir ${name}_microtrait  \\
        --nthreads ${task.cpus}
    """
}

// ============================================================
// PLASMIDFINDER  (Step 10a)
// Plasmid replicon detection via PlasmidFinder.
// Outputs JSON results; PLASMIDFINDER_TO_TSV will convert to TSV.
// Container: staphb/plasmidfinder  https://hub.docker.com/r/staphb/plasmidfinder
// Database bundled in container at /plasmidfinder_db by default;
// override with --plasmidfinder_db to use an external copy.
// ============================================================
process PLASMIDFINDER {
    tag "$name"
    publishDir "${params.outdir}/plasmidfinder_raw", mode: 'copy'
    container 'staphb/plasmidfinder:latest'
    conda     'bioconda::plasmidfinder'
    when: params.run_plasmidfinder

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}/results.json"), emit: json

    script:
    """
    mkdir -p ${name}
    python -m plasmidfinder \\
        -i ${fa_file}                      \\
        -o ${name}                         \\
        -p ${params.plasmidfinder_db}
    """
}

// ============================================================
// PLASMIDFINDER_TO_TSV
// Convert PlasmidFinder JSON output to TSV format.
// Schema: genome, contig, match_id, match_name, coverage, identity, hit_length
// ============================================================
process PLASMIDFINDER_TO_TSV {
    tag "$name"
    publishDir "${params.outdir}/plasmidfinder", mode: 'copy'
    container 'python:3.11-slim'
    conda     'conda-forge::python=3.11'
    when: params.run_plasmidfinder

    input:
    tuple val(name), path(results_json)

    output:
    tuple val(name), path("${name}_plasmidfinder.tsv")

    script:
    """
    python3 << 'EOF'
import json
import sys

# Load PlasmidFinder JSON
with open('${results_json}', 'r') as f:
    data = json.load(f)

# Write TSV with consistent schema
with open('${name}_plasmidfinder.tsv', 'w') as out:
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
// INTEGRONFINDER  (Step 10b)
// Detection of integrons and associated gene cassettes.
// Container: gem-pasteur/integron_finder  https://hub.docker.com/r/gempasteur/integron_finder
// ============================================================
process INTEGRONFINDER {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/integronfinder", mode: 'copy'
    container 'gempasteur/integron_finder:latest'
    conda     'bioconda::integron_finder'
    when: params.run_integronfinder

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}_integrons/"), emit: dir
    tuple val(name), path("${name}_integrons.tsv"), emit: tsv

    script:
    """
    integron_finder \\
        --outdir ${name}_integrons \\
        --cpu    ${task.cpus}      \\
        --local-max                \\
        ${fa_file}
    # Concatenate all .integrons tables into a single TSV (keep header from first file only)
    header=true
    for f in ${name}_integrons/Results_Integron_Finder_*/*.integrons; do
        if [ -f "\$f" ]; then
            if \$header; then
                cat "\$f" > ${name}_integrons.tsv
                header=false
            else
                grep -v '^#' "\$f" >> ${name}_integrons.tsv
            fi
        fi
    done
    touch ${name}_integrons.tsv
    """
}

// ============================================================
// MOB_SUITE  (Step 10c)
// Plasmid reconstruction, classification, and MOB typing.
// Container: kbessonov/mob_suite  https://hub.docker.com/r/kbessonov/mob_suite
// ============================================================
process MOB_SUITE {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/mob_suite", mode: 'copy'
    container 'kbessonov/mob_suite:latest'
    conda     'bioconda::mob_suite'
    when: params.run_mob_suite

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}_mob/"),   emit: dir
    tuple val(name), path("${name}_mob.tsv"), emit: tsv

    script:
    """
    mob_recon \\
        --infile      ${fa_file}   \\
        --outdir      ${name}_mob  \\
        --num_threads ${task.cpus} \\
        --run_typer
    cp ${name}_mob/contig_report.txt ${name}_mob.tsv
    """
}

// ============================================================
// ISESCAN  (Step 10d)
// Insertion sequence (IS) element detection and annotation.
// Container: staphb/isescan  https://hub.docker.com/r/staphb/isescan
// ============================================================
process ISESCAN {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/isescan", mode: 'copy'
    container 'staphb/isescan:latest'
    conda     'bioconda::isescan'
    when: params.run_isescan

    input:
    tuple val(name), path(fa_file)

    output:
    tuple val(name), path("${name}_isescan/"),   emit: dir
    tuple val(name), path("${name}_isescan.tsv"), emit: tsv

    script:
    """
    isescan.py \\
        --seqfile ${fa_file}          \\
        --output  ${name}_isescan     \\
        --nthread ${task.cpus}
    # Copy the main prediction table
    cp ${name}_isescan/*.tsv ${name}_isescan.tsv 2>/dev/null || touch ${name}_isescan.tsv
    """
}

// ============================================================
// ANTISMASH  (Step 11)
// Secondary metabolite biosynthetic gene cluster prediction.
// Container: antismash/standalone  https://hub.docker.com/r/antismash/standalone
// ============================================================
process ANTISMASH {
    tag "$name"
    label 'medium_cpu'
    publishDir "${params.outdir}/antismash", mode: 'copy'
    container 'antismash/standalone:latest'
    conda     'bioconda::antismash'
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
// Workflow
// ============================================================
workflow {

    // Collect, normalise, and (if needed) assemble all inputs into
    // (name, fa_file) per sample.  See subworkflows/input.nf for details.
    COLLECT_INPUTS()

    // Optionally filter by genome quality; if no report is given all samples proceed.
    def fa_ch
    if (params.quality_report) {
        def passing_names = Channel
            .fromPath(params.quality_report)
            .splitCsv(sep: '\t', skip: 1)
            .filter { row ->
                row[1].toFloat() >= params.min_completeness.toFloat() &&
                row[2].toFloat() <= params.max_contamination.toFloat()
            }
            .map { row -> [row[0], true] }
        fa_ch = COLLECT_INPUTS.out.samples
            .join(passing_names)
            .map { name, fa, _flag -> [name, fa] }
    } else {
        fa_ch = COLLECT_INPUTS.out.samples
    }

    def gtdbtk_db_ch
    if (params.run_taxonomy) {
        if (params.gtdbtk_db) {
            gtdbtk_db_ch = Channel.value(file(params.gtdbtk_db))
        } else {
            DOWNLOAD_GTDBTK_DB()
            gtdbtk_db_ch = DOWNLOAD_GTDBTK_DB.out.db.first()
        }
    }

    def bakta_db_ch
    if (params.run_bakta) {
        if (params.bakta_db) {
            bakta_db_ch = Channel.value(file(params.bakta_db))
        }
        else {
            DOWNLOAD_BAKTA_DB()
            bakta_db_ch = DOWNLOAD_BAKTA_DB.out.db.first()
        }
    }

    def kofam_db_ch
    if (params.run_kofam) {
        if (params.kofam_db) {
            kofam_db_ch = Channel.value(file(params.kofam_db))
        } else {
            DOWNLOAD_KOFAM_DB()
            kofam_db_ch = DOWNLOAD_KOFAM_DB.out.db.first()
        }
    }

    def eggnog_db_ch
    if (params.run_eggnog) {
        if (params.eggnog_db) {
            eggnog_db_ch = Channel.value(file(params.eggnog_db))
        } else {
            DOWNLOAD_EGGNOG_DB()
            eggnog_db_ch = DOWNLOAD_EGGNOG_DB.out.db.first()
        }
    }

    def card_db_ch
    if (params.run_card_rgi) {
        if (params.card_db) {
            card_db_ch = Channel.value(file(params.card_db))
        } else {
            DOWNLOAD_CARD_DB()
            card_db_ch = DOWNLOAD_CARD_DB.out.db.first()
        }
    }

    // Step 1 — Taxonomy: collect all FAs into one GTDB-Tk classify_wf run
    TAXONOMY(
        fa_ch.map { name, fa -> fa }.collect(),
        gtdbtk_db_ch
    )

    // Step 2 — Bakta: per-genome structural + functional annotation
    BAKTA(fa_ch, bakta_db_ch)

    // Steps 3, 5, 6, 8 — depend on Bakta protein sequences
    EGGNOG(BAKTA.out.faa, eggnog_db_ch)
    CARD_RGI(BAKTA.out.faa, card_db_ch)
    KOFAM(BAKTA.out.faa, kofam_db_ch)
    METABOLISMHMM(BAKTA.out.faa)

    // Step 6b — Reformat KOFAM output from space-separated to proper TSV
    KOFAM_REFORMAT(KOFAM.out.raw)

    // Steps 4, 7, 9, 10a-d, 11 — independent per-genome, run in parallel
    AMRFINDER(fa_ch)
    METABOLIC(fa_ch)
    MICROTRAIT(fa_ch)
    PLASMIDFINDER(fa_ch)
    INTEGRONFINDER(fa_ch)
    MOB_SUITE(fa_ch)
    ISESCAN(fa_ch)
    ANTISMASH(fa_ch)

    // Step 10a-b — Convert PlasmidFinder JSON to TSV
    PLASMIDFINDER_TO_TSV(PLASMIDFINDER.out.json)

    // Merge all TSVs into a Parquet file, then generate HTML report
    if (params.run_merge) {
        amrfinder_tsvs      = AMRFINDER.out.map            { name, tsv -> tsv }.collect()
        kofam_tsvs          = KOFAM_REFORMAT.out.map       { name, tsv -> tsv }.collect()
        card_rgi_tsvs       = CARD_RGI.out.tsv.map         { name, tsv -> tsv }.collect()
        plasmidfinder_tsvs  = PLASMIDFINDER_TO_TSV.out.map { name, tsv -> tsv }.collect()
        integronfinder_tsvs = INTEGRONFINDER.out.tsv.map   { name, tsv -> tsv }.collect()
        mob_suite_tsvs      = MOB_SUITE.out.tsv.map        { name, tsv -> tsv }.collect()
        isescan_tsvs        = ISESCAN.out.tsv.map          { name, tsv -> tsv }.collect()

        MERGE(
            amrfinder_tsvs,
            kofam_tsvs,
            card_rgi_tsvs,
            plasmidfinder_tsvs,
            integronfinder_tsvs,
            mob_suite_tsvs,
            isescan_tsvs
        )

        REPORT(MERGE.out.parquet)
    }
}
