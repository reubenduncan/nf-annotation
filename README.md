# nf-annotation

Nextflow pipeline for comprehensive annotation of metagenome-assembled genomes (MAGs). Accepts FASTA, FASTQ, or BAM inputs and runs up to 14 annotation tools in parallel, merging results into a single Parquet file and an HTML report.

## Pipeline overview

```
params.input/
    │
    ├─ FASTA (.fa/.fasta/.fna, plain or .gz)  → NORMALIZE_FASTA
    ├─ FASTQ (.fastq/.fq, plain or .gz)       → CONCAT_READS
    └─ BAM (.bam)                             → BAM_TO_FASTQ → CONCAT_READS
         │
         ├─[quality filter, optional]
         │
         ├──► TAXONOMY          (GTDB-Tk classify_wf)         → results/taxonomy/
         │
         ├──► BAKTA             (structural annotation)        → results/bakta/
         │       ├──► EGGNOG   (orthology/function)           → results/eggnog/
         │       ├──► CARD_RGI (AMR via CARD RGI)             → results/card_rgi/
         │       ├──► KOFAM    (KEGG orthology)
         │       │     └──► KOFAM_REFORMAT (space→TSV)        → results/kofam/
         │       └──► METABOLISHMM (HMM metabolic markers)    → results/metabolishmm/
         │
         ├──► AMRFINDER         (NCBI AMRFinderPlus)          → results/amrfinder/
         ├──► METABOLIC         (pathway reconstruction)       → results/metabolic/
         ├──► MICROTRAIT        (functional trait inference)   → results/microtrait/
         ├──► PLASMIDFINDER     (plasmid replicons)
         │     └──► PLASMIDFINDER_TO_TSV (JSON→TSV)           → results/plasmidfinder/
         ├──► INTEGRONFINDER    (integrons)                    → results/integronfinder/
         ├──► MOB_SUITE         (plasmid reconstruction/MOB)  → results/mob_suite/
         ├──► ISESCAN           (insertion sequences)         → results/isescan/
         └──► ANTISMASH         (biosynthetic gene clusters)  → results/antismash/
              │
              └──► MERGE        (all TSVs → Parquet)          → results/combined_results.parquet
                    └──► REPORT (HTML annotation report)      → results/annotation_report.html
```

BAKTA-dependent steps start as soon as each genome's Bakta job completes. All other per-genome steps run in parallel from the start.

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04
- Conda or Mamba (default runtime), **or** Docker / Singularity

## Runtime

The pipeline defaults to **conda** (mamba is used automatically when available). Override with a profile flag:

| Profile | Runtime |
|---------|---------|
| *(none)* | Conda / Mamba (default) |
| `-profile conda` | Conda / Mamba (explicit) |
| `-profile docker` | Docker |
| `-profile singularity` | Singularity |
| `-profile slurm,docker` | SLURM + Docker |
| `-profile slurm,singularity` | SLURM + Singularity |

## Input modes

The pipeline auto-detects the layout of `--input`:

| Layout | Behaviour |
|--------|-----------|
| **Flat** — files directly in `--input` | Each file (or R1/R2 pair) becomes one sample |
| **Subfolder** — immediate subdirectories in `--input` | All files in a subdir are merged into one sample |

Supported file types: `.fa`, `.fasta`, `.fna` (plain or `.gz`), `.fastq`, `.fq` (plain or `.gz`), `.bam`.

FASTQ and BAM inputs are converted to FASTA (reads treated as contigs — suitable for basecaller consensus reads). For short-read assemblies, assemble before running this pipeline.

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--input` | Path to input directory (FASTA / FASTQ / BAM files) |

### Databases

All database parameters are optional — if omitted the pipeline will auto-download the database to `--db_dir` (default: `databases/`) and cache it with `storeDir`.

| Parameter | Database | Notes |
|-----------|----------|-------|
| `--bakta_db` | [Bakta database](https://github.com/oschwengers/bakta#database) | ~30 GB |
| `--gtdbtk_db` | [GTDB-Tk reference data](https://ecogenomics.github.io/GTDBTk/) | ~66 GB |
| `--kofam_db` | [KofamScan profiles + ko_list](https://www.genome.jp/tools/kofamkoala/) | |
| `--eggnog_db` | [EggNOG-mapper data](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2) | |
| `--card_db` | [CARD database](https://card.mcmaster.ca/) | |
| `--plasmidfinder_db` | PlasmidFinder database | Defaults to `/plasmidfinder_db` bundled in container |

### Optional parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--quality_report` | `null` | TSV with columns: name, completeness, contamination. If provided, filters genomes before annotation |
| `--min_completeness` | `80` | Minimum completeness (%) to pass quality filter |
| `--max_contamination` | `5` | Maximum contamination (%) to pass quality filter |
| `--db_dir` | `databases` | Directory for auto-downloaded databases |
| `--outdir` | `results` | Output directory |

### Pipeline stages

All stages are enabled by default. Set any to `false` to skip:

| Parameter | Default | Tool |
|-----------|---------|------|
| `--run_taxonomy` | `true` | GTDB-Tk taxonomic classification |
| `--run_bakta` | `true` | Bakta structural annotation |
| `--run_eggnog` | `true` | EggNOG-mapper (requires `--run_bakta`) |
| `--run_amrfinder` | `true` | NCBI AMRFinderPlus |
| `--run_card_rgi` | `true` | CARD RGI AMR annotation (requires `--run_bakta`) |
| `--run_kofam` | `true` | KofamScan KEGG orthology (requires `--run_bakta`) |
| `--run_metabolic` | `true` | METABOLIC pathway reconstruction |
| `--run_metabolishmm` | `true` | metabolisHMM marker detection (requires `--run_bakta`) |
| `--run_microtrait` | `true` | microTrait functional trait inference |
| `--run_plasmidfinder` | `true` | PlasmidFinder plasmid detection |
| `--run_integronfinder` | `true` | IntegronFinder integron detection |
| `--run_mob_suite` | `true` | mob_suite plasmid reconstruction/MOB typing |
| `--run_isescan` | `true` | ISEScan insertion sequence detection |
| `--run_antismash` | `true` | antiSMASH biosynthetic gene cluster prediction |
| `--run_merge` | `true` | Merge all TSVs into Parquet |
| `--run_report` | `true` | Generate HTML annotation report |

## Usage

### Minimal example (pre-built databases)

```bash
nextflow run main.nf \
    --input          /path/to/genomes \
    --bakta_db       /path/to/bakta/db \
    --gtdbtk_db      /path/to/gtdbtk_data \
    --kofam_db       /path/to/kofam \
    --eggnog_db      /path/to/eggnog_data \
    --card_db        /path/to/card \
    --outdir         results
```

### Auto-download all databases

```bash
nextflow run main.nf \
    --input   /path/to/genomes \
    --db_dir  /path/to/databases \
    --outdir  results
```

### With quality filtering

```bash
nextflow run main.nf \
    --input              /path/to/genomes \
    --quality_report     /path/to/checkm_report.tsv \
    --min_completeness   80 \
    --max_contamination  5 \
    --outdir             results
```

### Skip specific stages

```bash
nextflow run main.nf \
    --input            /path/to/genomes \
    --run_metabolic    false \
    --run_metabolishmm false \
    --run_microtrait   false \
    --outdir           results
```

### Docker

```bash
nextflow run main.nf -profile docker [...]
```

### HPC / Singularity

```bash
nextflow run main.nf -profile slurm,singularity [...]
```

### HPC / conda

```bash
nextflow run main.nf -profile slurm [...]
```

### Resume a failed run

```bash
nextflow run main.nf -resume [...]
```

## Output structure

```
results/
├── taxonomy/              # GTDB-Tk classify_wf outputs
├── bakta/                 # Bakta per-genome outputs
├── eggnog/                # EggNOG-mapper TSVs
├── amrfinder/             # AMRFinderPlus TSVs
├── card_rgi/              # CARD RGI TSVs
├── kofam_raw/             # Raw KofamScan output (space-separated)
├── kofam/                 # Reformatted KofamScan TSVs
├── metabolic/             # METABOLIC per-genome outputs
├── metabolishmm/          # metabolisHMM per-genome outputs
├── microtrait/            # microTrait per-genome outputs
├── plasmidfinder_raw/     # Raw PlasmidFinder JSON
├── plasmidfinder/         # PlasmidFinder TSVs
├── integronfinder/        # IntegronFinder outputs + TSVs
├── mob_suite/             # mob_suite outputs + TSVs
├── isescan/               # ISEScan outputs + TSVs
├── antismash/             # antiSMASH per-genome outputs
├── combined_results.parquet   # Merged annotation results
└── annotation_report.html     # Self-contained HTML report
```

## Containers

| Step | Tool | Container |
|------|------|-----------|
| Input | NORMALIZE_FASTA | `ubuntu:24.04` |
| Input | BAM_TO_FASTQ | `staphb/samtools:latest` |
| Input | CONCAT_READS | `ubuntu:24.04` |
| 1 | GTDB-Tk | `nanozoo/gtdbtk:2.4.0--02c00d5` |
| 2 | Bakta | `oschwengers/bakta:latest` |
| 3 | EggNOG-mapper | `nanozoo/eggnog-mapper:2.1.13--c16a7d2` |
| 4 | AMRFinderPlus | `staphb/ncbi-amrfinderplus:latest` |
| 5 | CARD RGI | `finlaymaguire/rgi:latest` |
| 6 | KofamScan | `reubenduncan/kofam_scan:amd64` |
| 6b | KOFAM_REFORMAT | `python:3.11-slim` |
| 7 | METABOLIC | `jolespin/metabolic:v4.0` |
| 8 | metabolisHMM | `elizabethmcd/metabolishmm:latest` |
| 9 | microTrait | `ukaraoz/microtrait:latest` |
| 10a | PlasmidFinder | `staphb/plasmidfinder:latest` |
| 10a-b | PLASMIDFINDER_TO_TSV | `python:3.11-slim` |
| 10b | IntegronFinder | `gempasteur/integron_finder:latest` |
| 10c | mob_suite | `kbessonov/mob_suite:latest` |
| 10d | ISEScan | `staphb/isescan:latest` |
| 11 | antiSMASH | `antismash/standalone:latest` |
| Merge | MERGE | `python:3.11-slim` |
| Report | GENERATE_REPORT | `quay.io/biocontainers/pandas:1.5.2` |

## Data processing details

### Quality filtering

If `--quality_report` is provided, the pipeline reads the TSV (col1 = genome name, col2 = completeness, col3 = contamination), skipping one header row. Only genomes meeting `--min_completeness` and `--max_contamination` thresholds are passed to annotation tools.

### KOFAM reformatting

KofamScan outputs space-separated files (not true TSV). `KOFAM_REFORMAT` skips comment lines and writes proper tab-separated output:

```
gene_name    KO       threshold  score   E_value  KO_definition
APMBCM_00001 K21874   325.13     18.0    0.00031  SUN domain-containing protein 3
```

### PlasmidFinder JSON to TSV

`PLASMIDFINDER_TO_TSV` parses PlasmidFinder JSON and writes a TSV with consistent schema:

```
genome    contig       match_id   match_name         coverage  identity  hit_length
bin_001   contig_123   pDNA_001   plasmid_marker_A   95.2      98.5      1524
```

### Parquet merge

`MERGE` reads all tool-specific TSVs, adds a `source_tool` column, and concatenates them into a single Parquet file covering: `amrfinder`, `kofam`, `card_rgi`, `plasmidfinder`, `integronfinder`, `mob_suite`, `isescan`.
