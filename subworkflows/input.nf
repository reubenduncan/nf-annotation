// subworkflows/input.nf
//
// Reusable input-collection subworkflow.
//
// Scans params.input for sequence files, groups them into per-sample channels,
// and emits sanitized (name, fa_file) tuples ready for downstream annotation.
//
// Supported file types:
//   FASTA  : .fa  .fasta  .fna  (plain or .gz)
//   FASTQ  : .fastq  .fq        (plain or .gz)  → assembled with MEGAHIT
//   BAM    : .bam               → converted to FASTQ, then assembled
//
// Input modes (auto-detected from params.input layout):
//   Flat mode      — files sit directly in params.input
//                    each file (or R1/R2 pair) becomes one sample
//   Subfolder mode — immediate subdirectories of params.input
//                    all files inside a subdir are merged into one sample

// ============================================================
// NORMALIZE_FASTA
// Decompresses and concatenates any mix of FASTA files (.fa/.fasta/.fna,
// plain or .gz) into a single uncompressed .fa for a sample.
// ============================================================
process NORMALIZE_FASTA {
    tag "$name"
    container 'ubuntu:24.04'

    input:
    tuple val(name), path(fasta_files)

    output:
    tuple val(name), path("${name}.fa")

    script:
    """
    for f in ${fasta_files}; do
        case "\$f" in
            *.gz) gunzip -c "\$f" ;;
            *)    cat "\$f" ;;
        esac
    done > ${name}.fa
    """
}

// ============================================================
// BAM_TO_FASTQ
// Merges multiple BAMs if needed, name-sorts, and converts to
// paired FASTQ.  Empty output files (e.g. no unpaired reads)
// are removed so ASSEMBLE_READS sees only non-empty inputs.
// ============================================================
process BAM_TO_FASTQ {
    tag "$name"
    label 'medium_cpu'
    container 'staphb/samtools:latest'

    input:
    tuple val(name), path(bam_files)

    output:
    tuple val(name), path("reads_*.fastq.gz")

    script:
    """
    samtools cat ${bam_files} \
        | samtools sort -n -@ ${task.cpus} \
        | samtools fastq -@ ${task.cpus} \
            -1 reads_R1.fastq.gz \
            -2 reads_R2.fastq.gz \
            -s reads_unpaired.fastq.gz \
            -0 /dev/null \
            -
    find . -name 'reads_*.fastq.gz' -empty -delete
    """
}

// ============================================================
// CONCAT_READS
// Concatenates FASTQ files and converts to FASTA, treating each
// read as a contig directly (no assembly).  Handles both plain
// and gzipped input.  Suitable for basecaller consensus reads
// where each read is already a complete sequence.
// ============================================================
process CONCAT_READS {
    tag "$name"
    container 'ubuntu:24.04'

    input:
    tuple val(name), path(fastq_files)

    output:
    tuple val(name), path("${name}.fa")

    script:
    """
    for f in ${fastq_files}; do
        case "\$f" in
            *.gz) gunzip -c "\$f" ;;
            *)    cat "\$f" ;;
        esac
    done | awk 'NR%4==1{print ">"substr(\$0,2)} NR%4==2{print}' > ${name}.fa
    """
}

// ============================================================
// COLLECT_INPUTS
//
// Emits:
//   samples — channel of (name, fa_file) tuples
// ============================================================
workflow COLLECT_INPUTS {

    main:
    def dir        = file(params.input)
    def has_subdirs = dir.listFiles()?.any { it.isDirectory() } ?: false

    def all_files_ch
    if (has_subdirs) {
        // Subfolder mode: each immediate subdirectory = one sample
        all_files_ch = Channel
            .fromPath("${params.input}/*/*")
            .filter { f ->
                f.name =~ /\.(fa|fasta|fna)(\.gz)?$/ ||
                f.name =~ /\.(fastq|fq)(\.gz)?$/     ||
                f.name.endsWith('.bam')
            }
            .map { f -> [f.parent.name, f] }
    } else {
        // Flat mode: one sample per file.
        // R1/R2 FASTQ pairs are grouped by stripping pairing suffixes.
        all_files_ch = Channel
            .fromPath("${params.input}/*")
            .filter { f ->
                f.name =~ /\.(fa|fasta|fna)(\.gz)?$/ ||
                f.name =~ /\.(fastq|fq)(\.gz)?$/     ||
                f.name.endsWith('.bam')
            }
            .map { f ->
                def stem = f.name
                    .replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
                    .replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                    .replaceAll(/\.bam$/, '')
                    .replaceAll(/_R[12](_\d+)?$/, '')
                    .replaceAll(/_[12]$/, '')
                [stem, f]
            }
    }

    // Group files by sample, classify by dominant type, route to processes.
    // Priority when a sample has mixed types: bam > fastq > fasta.
    def classified = all_files_ch
        .groupTuple()
        .map { name, files ->
            def type = files.any { it.name.endsWith('.bam') }           ? 'bam'
                     : files.any { it.name =~ /\.(fastq|fq)(\.gz)?$/ } ? 'fastq'
                     :                                                    'fasta'
            [name, type, files]
        }
        .branch {
            fasta: it[1] == 'fasta'
            fastq: it[1] == 'fastq'
            bam:   it[1] == 'bam'
        }

    NORMALIZE_FASTA(classified.fasta.map { name, type, files -> [name, files] })
    BAM_TO_FASTQ   (classified.bam.map   { name, type, files -> [name, files] })

    // Concatenate + convert both direct FASTQ inputs and BAM-derived FASTQ
    CONCAT_READS(
        classified.fastq.map { name, type, files -> [name, files] }
            .mix(BAM_TO_FASTQ.out)
    )

    emit:
    samples = NORMALIZE_FASTA.out.mix(CONCAT_READS.out)
}
