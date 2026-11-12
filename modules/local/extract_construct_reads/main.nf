// Step 2b: cut the two small BAMs out of the full alignment.
//
//   <id>.construct_reads.bam  every alignment of every read touching the construct
//   <id>.split.bam            of those, the reads broken across two places
//
// Read NAMES are collected on the construct first, then ALL alignments of those
// reads are pulled: that brings the genome-side supplementary alignments with
// them. Asking for the construct region directly would return only the
// construct half of each read.
//
// These are the files to open in IGV, and every analysis below reads them
// rather than the full BAM.

process EXTRACT_CONSTRUCT_READS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.construct_reads.bam"), path("${meta.id}.construct_reads.bam.bai"), emit: construct_reads
    tuple val(meta), path("${meta.id}.split.bam"), path("${meta.id}.split.bam.bai"),                     emit: split
    tuple val(meta), path("${meta.id}.names.txt"),                                                       emit: names
    tuple val(meta), path("${meta.id}.contigs_hit.txt"),                                                 emit: contigs
    tuple val(meta), path("${meta.id}.construct_reads.sam.gz"),                                          emit: sam

    script:
    def construct = meta.construct_contig
    """
    set -euo pipefail

    samtools view ${bam} '${construct}' | cut -f1 | sort -u > ${meta.id}.names.txt

    samtools view -b -N ${meta.id}.names.txt -o ${meta.id}.construct_reads.bam ${bam}
    samtools index -@ ${task.cpus} ${meta.id}.construct_reads.bam

    samtools view -b -d SA -o ${meta.id}.split.bam ${meta.id}.construct_reads.bam
    samtools index -@ ${task.cpus} ${meta.id}.split.bam

    # Kept for the optional architecture step, which runs in a python container
    # with no samtools in it.
    samtools view -F 0x104 ${meta.id}.construct_reads.bam | gzip > ${meta.id}.construct_reads.sam.gz

    {
        echo "construct primary reads: \$(samtools view -c -F 0x900 ${bam} '${construct}')"
        echo "total primary reads:     \$(samtools view -c -F 0x904 ${bam})"
        echo "reads touching construct: \$(wc -l < ${meta.id}.names.txt | tr -d ' ')"
        echo
        echo "--- non-construct contigs hit by these reads ---"
        # -F 0x104 drops unmapped + secondary but KEEPS supplementary (0x800).
        # On a junction read the construct is usually the primary and the genome
        # side the supplementary, so -F 0x904 would hide half the evidence.
        samtools view -F 0x104 ${meta.id}.split.bam \\
            | awk -v construct='${construct}' '\$3 != construct {print \$3}' \\
            | sort | uniq -c | sort -rn
    } > ${meta.id}.contigs_hit.txt

    cat ${meta.id}.contigs_hit.txt
    """

    stub:
    """
    touch ${meta.id}.construct_reads.bam ${meta.id}.construct_reads.bam.bai
    touch ${meta.id}.split.bam ${meta.id}.split.bam.bai
    touch ${meta.id}.names.txt ${meta.id}.contigs_hit.txt
    echo | gzip > ${meta.id}.construct_reads.sam.gz
    """
}
