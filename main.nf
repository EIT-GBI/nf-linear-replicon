#!/usr/bin/env nextflow

// nf-linear-replicon - where a linear replicon construct integrated into the
// host genome, from ONT reads.
//
// The whole idea in one line: align each sample against genome AND construct
// TOGETHER, so a read crossing an integration junction gets a primary
// alignment on one sequence and a supplementary on the other. Every result
// below is read off those pairs.

//
// Modules loaded from nf-mod-* submodules
//
include { DORADO_ALIGNER     } from './modules/dorado/aligner/main.nf'
include { MINIMAP2_ALIGN     } from './modules/minimap2/align/main.nf'
include { SAMTOOLS_FAIDX     } from './modules/samtools/faidx/main.nf'
include { SAMTOOLS_FLAGSTAT  } from './modules/samtools/flagstat/main.nf'

//
// Modules local to this pipeline
//
include { PREPARE_READS            } from './modules/local/prepare_reads/main.nf'
include { BUILD_REFERENCE          } from './modules/local/build_reference/main.nf'
include { CONSTRUCT_HOMOLOGY       } from './modules/local/construct_homology/main.nf'
include { EXTRACT_CONSTRUCT_READS  } from './modules/local/extract_construct_reads/main.nf'
include { JUNCTION_POSITIONS       } from './modules/local/junction_positions/main.nf'
include { JUNCTION_REPORT          } from './modules/local/junction_report/main.nf'
include { JUNCTION_MATRIX          } from './modules/local/junction_matrix/main.nf'
include { READ_ARCHITECTURE        } from './modules/local/read_architecture/main.nf'


// A reference named in the samplesheet, resolved under params.reference_dir.
// The name may carry its extension or not: the cluster samplesheets write
// bare names ('CHOWT_genome'), so the same extensions are tried here.
def resolveReference(String name) {
    def given = file("${params.reference_dir}/${name}")
    if (given.exists())
        return given

    def extensions = ['.fa', '.fasta', '.fna', '.fa.gz', '.fasta.gz', '.fna.gz']
    def found = extensions
        .collect { ext -> file("${params.reference_dir}/${name}${ext}") }
        .find { candidate -> candidate.exists() }
    if (found)
        return found

    error "No FASTA for '${name}' in ${params.reference_dir} (tried the name as given, and .fa/.fasta/.fna, plain and .gz)"
}

// Name of the genome+construct pair, used to key the reference channels and to
// name the built files. Basenames only: the samplesheet may give sub-paths.
def pairTag(genome, construct) {
    def genome_name    = genome.name.replaceFirst(/\.(fa|fasta|fna)(\.gz)?$/, '')
    def construct_name = construct.name.replaceFirst(/\.(fa|fasta|fna)(\.gz)?$/, '')
    return "${genome_name}__${construct_name}"
}

// One samplesheet column, matched the way the cluster scripts matched it:
// case-insensitive, trimmed, CRLF- and BOM-safe, and BY PREFIX -- so
// `additional_reference`, `additional_reference_name` and
// `additional_reference_sequence` are all the same column. Instrument-written
// samplesheets vary, and an unmatched column reaching the pipeline as null is
// how this went wrong before.
def normaliseHeader(key) {
    return key?.toString()?.replace('\uFEFF', '')?.trim()?.toLowerCase()?.replaceAll(/[^a-z0-9]/, '') ?: ''
}

def column(Map row, String want) {
    def wanted = want.toLowerCase().replaceAll(/[^a-z0-9]/, '')

    def hit = row.find { key, _value -> normaliseHeader(key) == wanted }
    if (!hit)
        hit = row.find { key, _value -> normaliseHeader(key).startsWith(wanted) }

    def value = hit?.value?.toString()?.replace('\r', '')?.trim()
    return value ?: null
}

// The same, but a missing column is a hard error naming the headers that are
// there -- far easier to act on than a null three steps downstream.
def requiredColumn(Map row, String want) {
    def value = column(row, want)
    if (!value)
        error "Samplesheet has no '${want}' column, or the value is empty in one row.\n" +
              "       Columns found: ${row.keySet().join(', ')}"
    return value
}

// Sample names come from instrument samplesheets and carry spaces, brackets,
// '+' and '%'. They end up in file names and on command lines, so they are
// reduced to a safe token here. The original is kept in meta.sample_name.
def safeName(String name) {
    def safe = name.replaceAll(/[^A-Za-z0-9._-]+/, '_').replaceAll(/^_+/, '').replaceAll(/_+$/, '')
    return safe ?: 'unnamed'
}

// The reads for one samplesheet row. fastq_path is a FOLDER, and may hold
// either FASTQs (fastq_pass) or unaligned BAMs (bam_pass, what dorado writes
// by default now). FASTQ wins if a folder somehow holds both.
//
// A relative fastq_path is taken relative to the samplesheet itself, so a
// samplesheet and the reads it names can be moved together.
def readsIn(String path, base) {
    // The test is on the STRING: file() resolves a relative path against the
    // launch directory, so the Path it hands back is always absolute.
    def dir = (path.startsWith('/') || path.contains('://')) ? file(path) : file("${base}/${path}")
    if (!dir.exists())
        error "fastq_path does not exist: ${dir}"

    // files(), not file(): file() with a glob returns a bare Path when exactly
    // one file matches, and a Path iterates over its own name components.
    def fastqs = (files("${dir}/**.fastq.gz") + files("${dir}/**.fastq")).sort { fastq -> fastq.name }
    if (fastqs)
        return [type: 'fastq', files: fastqs]

    def bams = files("${dir}/**.bam").sort { bam -> bam.name }
    if (bams)
        return [type: 'bam', files: bams]

    error "No *.fastq.gz, *.fastq or *.bam under ${dir}"
}


workflow {

    main:

    if (!params.reference_dir)
        error "You must set params.reference_dir: the folder holding the genome and construct FASTAs."
    if (!params.samplesheet)
        error "You must set params.samplesheet."

    // A folder here is the easy mistake -- the run folder rather than the CSV
    // inside it -- and without this check it surfaces as a Java stack trace
    // from splitCsv ("Is a directory") with nothing naming the parameter.
    def samplesheet_path = file(params.samplesheet)
    if (samplesheet_path.isDirectory())
        error "params.samplesheet is a folder, not a CSV file: ${params.samplesheet}\n" +
              "       Point it at the samplesheet itself, e.g. ${params.samplesheet}/samplesheet.csv"
    if (!(params.alignment.tool in ['minimap2', 'dorado']))
        error "Invalid value for params.alignment.tool: ${params.alignment.tool}. Use 'minimap2' or 'dorado'."

    samplesheet_ch = channel.fromPath(params.samplesheet, checkIfExists: true)
    def samplesheet_dir = file(params.samplesheet).parent

    // One entry per samplesheet row. Everything the run needs from the row is
    // resolved here, so a bad path fails in seconds rather than hours in.
    rows_ch = samplesheet_ch
        .splitCsv(header: true, strip: true)
        // drop entirely blank rows (trailing ",,," lines are common in exported CSVs)
        .filter { row -> row.values().any { v -> v?.toString()?.trim() } }
        .map { row ->
            def sample    = requiredColumn(row, 'sample')
            def genome    = resolveReference(requiredColumn(row, 'reference_genome'))
            def construct = resolveReference(requiredColumn(row, 'additional'))

            // Sample names repeat across ReadSets (same clone, two runs), so
            // outputs are keyed on all three or they collide.
            def readset = column(row, 'readset') ?: 'readset'
            def barcode = column(row, 'barcode') ?: 'barcode'

            def meta = [
                id          : "${safeName(sample)}__${safeName(readset)}__${safeName(barcode)}",
                sample      : safeName(sample),
                sample_name : sample,
                readset     : readset,
                barcode     : barcode,
                tag         : pairTag(genome, construct),
                platform    : column(row, 'platform') ?: 'ONT',
                preset      : params.preset,
            ]
            tuple(meta, genome, construct, readsIn(requiredColumn(row, 'fastq_path'), samplesheet_dir))
        }

    // ---------------------------------------------------------------------
    // Step 1: one combined genome+construct reference per pair
    // ---------------------------------------------------------------------
    pairs_ch = rows_ch
        .map { meta, genome, construct, _reads -> tuple(meta.tag, genome, construct) }
        .unique()

    BUILD_REFERENCE(pairs_ch)
    SAMTOOLS_FAIDX(BUILD_REFERENCE.out.fasta.map { tag, fasta -> tuple([id: tag], fasta) })

    // Shared sequence between construct and genome produces fake junctions at a
    // FIXED construct coordinate in every sample. Read this report before
    // believing any result.
    CONSTRUCT_HOMOLOGY(pairs_ch)

    // tuple(tag, fasta, fai, construct contig name)
    reference_ch = BUILD_REFERENCE.out.fasta
        .join(BUILD_REFERENCE.out.contig)
        .join(SAMTOOLS_FAIDX.out.fai.map { meta, fai -> tuple(meta.id, fai) })
        .map { tag, fasta, contig, fai -> tuple(tag, fasta, fai, contig.text.trim()) }

    // ---------------------------------------------------------------------
    // Step 2: reads -> one sorted BAM per row
    // ---------------------------------------------------------------------
    // Every row's read files are checked and joined into ONE file first: both
    // aligners take a single reads argument, and a bam_pass/barcodeNN folder
    // holds many. For dorado that join is `samtools cat` (BGZF block copy, tags
    // intact); for minimap2 the BAMs have to be decoded to FASTQ.
    PREPARE_READS(rows_ch.map { meta, _genome, _construct, reads ->
        tuple(meta + [read_type: reads.type], reads.files)
    })

    align_in = PREPARE_READS.out.reads
        .map { meta, reads -> tuple(meta.tag, meta, reads) }
        .combine(reference_ch, by: 0)
        .multiMap { _tag, meta, reads, fasta, fai, construct ->
            reads: tuple(meta + [construct_contig: construct], reads)
            fasta: tuple(fasta, fai)
        }

    // Both aligners emit tuple(meta, bam, bai).
    if (params.alignment.tool == 'dorado') {
        DORADO_ALIGNER(align_in.reads, align_in.fasta)
        bam_ch = DORADO_ALIGNER.out.bam
    }
    else {
        // -Y keeps the split pieces soft-clipped, so the sequence either side
        // of a junction stays in the record; see ext.args in nextflow.config.
        MINIMAP2_ALIGN(align_in.reads, align_in.fasta)
        bam_ch = MINIMAP2_ALIGN.out.bam
    }

    SAMTOOLS_FLAGSTAT(bam_ch.map { meta, bam, _bai -> tuple(meta, bam) })

    // Every alignment of every read that touches the construct, and of those,
    // the reads broken across two places. These are the small BAMs to open in
    // IGV, and everything below reads them rather than the full BAM.
    EXTRACT_CONSTRUCT_READS(bam_ch)

    // ---------------------------------------------------------------------
    // Step 3a: per-sample junction positions, then one report across samples
    // ---------------------------------------------------------------------
    JUNCTION_POSITIONS(EXTRACT_CONSTRUCT_READS.out.split)
    JUNCTION_REPORT(JUNCTION_POSITIONS.out.csv.map { _meta, csv -> csv }.collect())

    // ---------------------------------------------------------------------
    // Step 3b: the same junctions pooled across samples. A real integration
    // recurs at one locus in independent libraries of the same clone; ONT
    // chimeras scatter and do not.
    // ---------------------------------------------------------------------
    matrix_in = EXTRACT_CONSTRUCT_READS.out.split
        .join(EXTRACT_CONSTRUCT_READS.out.names)
        .multiMap { _meta, split_bam, _bai, names ->
            split: split_bam
            names: names
        }
    JUNCTION_MATRIX(matrix_in.split.collect(), matrix_in.names.collect())

    // ---------------------------------------------------------------------
    // Optional: per-read architecture of the construct reads
    // ---------------------------------------------------------------------
    ch_architecture = channel.empty()
    if (params.architecture.enabled) {
        architecture_in = EXTRACT_CONSTRUCT_READS.out.sam
            .map { meta, sam -> tuple(meta.tag, meta, sam) }
            .combine(BUILD_REFERENCE.out.construct_fasta, by: 0)
            .map { _tag, meta, sam, construct_fasta -> tuple(meta, sam, construct_fasta) }

        READ_ARCHITECTURE(architecture_in, file("${projectDir}/" + params.scripts.py_architecture, checkIfExists: true))
        ch_architecture = READ_ARCHITECTURE.out.report
    }

    publish:
    samplesheet     = samplesheet_ch
    reference       = BUILD_REFERENCE.out.fasta.map { _tag, fasta -> fasta }
    reference_fai   = SAMTOOLS_FAIDX.out.fai
    homology        = CONSTRUCT_HOMOLOGY.out.report
    // DORADO_ALIGNER carries a publishDir of its own (it is a submodule, so
    // not ours to change) and already writes to <outdir>/aligned. Publishing
    // the same files from here as well would write each BAM twice.
    alignment       = params.alignment.tool == 'dorado' ? channel.empty() : bam_ch
    flagstat        = SAMTOOLS_FLAGSTAT.out.flagstat
    subset          = EXTRACT_CONSTRUCT_READS.out.construct_reads
    split           = EXTRACT_CONSTRUCT_READS.out.split
    names           = EXTRACT_CONSTRUCT_READS.out.names
    contigs_hit     = EXTRACT_CONSTRUCT_READS.out.contigs
    junctions       = JUNCTION_REPORT.out.csv
    distribution    = JUNCTION_REPORT.out.distribution
    candidates      = JUNCTION_REPORT.out.candidates
    matrix          = JUNCTION_MATRIX.out.csv
    matrix_summary  = JUNCTION_MATRIX.out.summary
    architecture    = ch_architecture
}

//
// Where each published output lands under the output directory
// (params.outdir, overridable with -output-dir / -o).
//
// Channels carrying tuple(meta, files...) are scanned recursively, so the meta
// map is ignored and every file in the tuple is published to the given path.
//
output {
    samplesheet   { path 'samplesheet'   }
    reference     { path 'ref'           }
    reference_fai { path 'ref'           }
    homology      { path 'ref'           }
    alignment     { path 'aligned'       }
    flagstat      { path 'qc/flagstat'   }
    subset        { path 'subset'        }
    split         { path 'subset'        }
    names         { path 'subset'        }
    contigs_hit   { path 'qc/contigs'    }
    junctions     { path 'results'       }
    distribution  { path 'results'       }
    candidates    { path 'results'       }
    matrix        { path 'results'       }
    matrix_summary { path 'results'      }
    architecture  { path 'results/architecture' }
}
