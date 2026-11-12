// One checked, single file of reads per samplesheet row.
//
// Two jobs, done in one pass because they are one read of the same files:
//
// 1. INTEGRITY. A truncated .fastq.gz desynchronises the SAM stream and
//    samtools aborts hours later with the misleading "CIGAR and query sequence
//    are of different length". `samtools quickcheck` is the BAM equivalent and
//    a stricter one: it reads each header and verifies the EOF block.
//      -u  bam_pass is UNALIGNED and has no @SQ lines; without this every file
//          is reported "had no targets in header" and the job dies before it
//          starts.
//
// 2. ONE FILE. `dorado aligner` takes a single reads argument. One barcode
//    folder often holds just one BAM, and then this step only renames it --
//    a symlink, no copy, no data read. When MinKNOW has written the barcode in
//    batches, `samtools cat` joins them by copying BGZF blocks: no decode, no
//    re-encode, and the MM/ML modified-base tags ride along untouched.
//
// For the minimap2 route the same files have to become FASTQ instead, because
// minimap2 cannot read BAM. That is the one case where the reads are decoded,
// and -T MM,ML (plus -y on the aligner) is what carries the tags across.

process PREPARE_READS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.reads.*"), emit: reads

    script:
    // dorado reads BAM directly; minimap2 needs FASTQ
    def target = params.alignment.tool == 'dorado' ? 'bam' : 'fastq'
    def args   = task.ext.args ?: '-F 0x900 -T MM,ML'
    """
    set -euo pipefail

    for f in ${reads}; do
        case "\$f" in
            *.bam) samtools quickcheck -v -u "\$f" || { echo "CORRUPT bam: \$f"   >&2; exit 1; } ;;
            *.gz)  gzip -t "\$f"                  || { echo "CORRUPT fastq: \$f" >&2; exit 1; } ;;
        esac
    done

    set -- ${reads}

    if [ "${meta.read_type}" = bam ] && [ "${target}" = bam ]; then
        if [ "\$#" -eq 1 ]; then
            ln -s "\$1" ${meta.id}.reads.bam      # nothing to join
        else
            samtools cat -o ${meta.id}.reads.bam ${reads}
        fi
    elif [ "${meta.read_type}" = bam ]; then
        # set -o pipefail makes a bad file fail the job loudly rather than
        # truncate the output silently. gzip -1: this file is read once.
        for bam in ${reads}; do
            samtools fastq ${args} -@ ${task.cpus} "\$bam"
        done | gzip -1 > ${meta.id}.reads.fastq.gz
    elif [ "\$#" -eq 1 ] && [ "\${1##*.}" = gz ]; then
        ln -s "\$1" ${meta.id}.reads.fastq.gz    # nothing to join
    else
        # gzip members concatenate, so joining .gz files is just cat
        for fastq in ${reads}; do
            case "\$fastq" in
                *.gz) cat "\$fastq" ;;
                *)    gzip -1 -c "\$fastq" ;;
            esac
        done > ${meta.id}.reads.fastq.gz
    fi

    echo "${meta.id}: \$(ls ${reads} | wc -l | tr -d ' ') ${meta.read_type} file(s) checked -> \$(ls ${meta.id}.reads.*)"
    """

    stub:
    """
    if [ "${params.alignment.tool}" = dorado ] && [ "${meta.read_type}" = bam ]; then
        touch ${meta.id}.reads.bam
    else
        echo | gzip > ${meta.id}.reads.fastq.gz
    fi
    """
}
