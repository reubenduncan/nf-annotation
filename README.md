# Gene Assembly Annotation Pipeline (Nextflow)

8-step metagenome-assembled genome (MAG) annotation pipeline, with optional reformatting and merging of all results into a single Parquet file.

## Pipeline overview

```
bins/*.fa.gz
    │
    ├─[quality filter]──► TAXONOMY    (GTDB-Tk)          → results/taxonomy/
    │
    └─[all bins]
         │
         ├──► BAKTA        (structural annotation)       → results/bakta/
         │       ├──► EGGNOG   (orthology/function)      → results/eggnog/
         │       └──► KOFAM    (KEGG orthology)
         │            └──► KOFAM_REFORMAT (space→TSV)   → results/kofam/
         │
         ├──► AMRFINDER    (AMR genes)                   → results/amrfinder/
         ├──► VFDB         (virulence factors)           → results/vfdb/
         ├──► MGE          (plasmids/mobile elements)
         │     └──► MGE_TO_TSV (JSON→TSV)               → results/mge/
         ├──► ANTISMASH    (biosynthetic gene clusters)  → results/antismash/
         │
         └──► MERGE_RESULTS (all TSVs→Parquet)          → results/combined_results.parquet
```

Steps 4–8 run in parallel. Steps 3 and 6 start as soon as each genome's Bakta job finishes. Reformatting steps (KOFAM_REFORMAT, MGE_TO_TSV) and final merge can be enabled/disabled independently.

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04
- Docker (or Singularity — see profiles below)
- Python with pandas library (for merge step, automatically pulled in `python:3.11-slim` container)

## Databases

Download these separately and provide their paths at runtime:

| Parameter | Database |
|-----------|----------|
| `--gtdbtk_db` | [GTDB-Tk reference data](https://ecogenomics.github.io/GTDBTk/installing/index.html) |
| `--bakta_db` | [Bakta database](https://github.com/oschwengers/bakta#database) |
| `--eggnog_db` | [EggNOG-mapper data](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2) (`download_eggnog_data.py`) |
| `--kofam_db` | [KofamScan profiles + ko_list](https://www.genome.jp/tools/kofamkoala/) |
| `--plasmidfinder_db` | Optional — defaults to the database bundled in the `staphb/plasmidfinder` container |

## Containers used

| Step | Tool | Container |
|------|------|-----------|
| 1 | GTDB-Tk | `ecogenomics/gtdbtk:2.4.0` |
| 2 | Bakta | `oschwengers/bakta:latest` |
| 3 | EggNOG-mapper | `nanozoo/eggnog-mapper:2.1.12--0` |
| 4 | AMRFinderPlus | `staphb/ncbi-amrfinderplus:latest` |
| 5 | abricate (VFDB) | `staphb/abricate:latest` |
| 6 | KofamScan | `reubenduncan/kofam_scan:amd64` |
| 6b | KOFAM_REFORMAT | `python:3.11-slim` |
| 7 | PlasmidFinder | `staphb/plasmidfinder:latest` |
| 7b | MGE_TO_TSV | `python:3.11-slim` |
| 8 | antiSMASH | `antismash/antismash:latest` |
| Merge | MERGE_RESULTS | `python:3.11-slim` |

## Usage

### Minimal example

```bash
nextflow run main.nf \
    --bins_dir       /path/to/final_bins \
    --quality_report /path/to/quality_report.tsv \
    --bakta_db       /path/to/bakta/db \
    --gtdbtk_db      /path/to/gtdbtk_data \
    --kofam_db       /path/to/kofam \
    --eggnog_db      /path/to/eggnog_data \
    --run_bakta \
    --run_amrfinder \
    --outdir         results
```

### Enable all stages including merge

```bash
nextflow run main.nf \
    --bins_dir       /path/to/final_bins \
    --quality_report /path/to/quality_report.tsv \
    --bakta_db       /path/to/bakta/db \
    --gtdbtk_db      /path/to/gtdbtk_data \
    --kofam_db       /path/to/kofam \
    --eggnog_db      /path/to/eggnog_data \
    --run_taxonomy \
    --run_bakta \
    --run_eggnog \
    --run_amrfinder \
    --run_vfdb \
    --run_kofam \
    --run_mge \
    --run_antismash \
    --run_merge \
    --outdir         results
```

### Skippable stages

Each annotation stage is disabled by default and can be enabled with `--run_<stage>`:

- `--run_taxonomy` → runs GTDB-Tk taxonomic classification
- `--run_bakta` → runs Bakta structural annotation
- `--run_eggnog` → runs EggNOG-mapper (requires `--run_bakta`)
- `--run_amrfinder` → runs AMRFinderPlus
- `--run_vfdb` → runs abricate with VFDB
- `--run_kofam` → runs KofamScan + KOFAM_REFORMAT (requires `--run_bakta`)
- `--run_mge` → runs PlasmidFinder + MGE_TO_TSV
- `--run_antismash` → runs antiSMASH
- `--run_merge` → merges all TSVs into a single Parquet file

### Quality filtering thresholds (optional)

Defaults match the original scripts (≥80% completeness, ≤5% contamination):

```bash
    --min_completeness 80 \
    --max_contamination 5
```

### HPC / Singularity

```bash
nextflow run main.nf -profile slurm,singularity [...]
```

### Resume a failed run

```bash
nextflow run main.nf -resume [...]
```

## Output structure

```
results/
├── taxonomy/              # GTDB-Tk outputs (if --run_taxonomy)
├── bakta/                 # Bakta per-genome outputs (if --run_bakta)
├── eggnog/                # EggNOG-mapper TSVs (if --run_eggnog)
├── amrfinder/             # AMRFinderPlus TSVs (if --run_amrfinder)
├── vfdb/                  # VFDB TSVs (if --run_vfdb)
├── kofam_raw/             # Raw space-separated KOFAM output
├── kofam/                 # Reformatted KOFAM TSVs (if --run_kofam)
├── mge_raw/               # Raw PlasmidFinder JSON
├── mge/                   # Reformatted MGE TSVs (if --run_mge)
├── antismash/             # antiSMASH per-genome outputs (if --run_antismash)
├── combined_results.parquet   # Merged results (if --run_merge)
├── pipeline_trace.txt     # Execution trace
└── pipeline_report.html   # HTML execution report
```

## Data processing details

### KOFAM reformatting
KofamScan outputs space-separated files (not true TSV). KOFAM_REFORMAT reads the raw output, skips comment lines, and writes proper tab-separated output with a header:

```
gene_name    KO       threshold  score   E_value  KO_definition
APMBCM_00001 K21874   325.13     18.0    0.00031  SUN domain-containing protein 3
```

### MGE JSON to TSV
PlasmidFinder outputs JSON. MGE_TO_TSV parses the results and writes a TSV with consistent schema:

```
genome    contig       match_id   match_name         coverage  identity  hit_length
bin_001   contig_123   pDNA_001   plasmid_marker_A   95.2      98.5      1524
```

### Parquet merge
MERGE_RESULTS reads all tool-specific TSVs, adds a `source_tool` column (one of: `amrfinder`, `vfdb`, `kofam`, `mge_plasmidfinder`), and concatenates them into a single Parquet file. It:
- Skips comment lines (lines starting with `#`)
- Handles `##` comment blocks gracefully
- Aligns columns by name, so different tools' outputs are preserved as-is with tool annotation
- Outputs to `results/combined_results.parquet`

## Quality report format

The TSV must have genome name in column 1, completeness in column 2, contamination in column 3 (no header, or with one header row — the pipeline skips row 1 automatically):

```
genome_name    completeness    contamination
bin_001        95.3            1.2
bin_002        72.1            8.4
```


