// Optional: the architecture of each construct read, and the integration loci
// those reads define.
//
// Every read is one independent observation of one molecule, so nothing aligns
// reads to each other -- the structure is already in the alignments. A real
// integration site is FIXED on the junction side and SCATTERED on the shear
// side; distinct shear points count independent molecules, which duplicates and
// duplex read-pairs cannot fake.
//
// architecture.py is taken verbatim from the cluster pipeline. It reads
// `samtools view -F 0x104` output on stdin, which EXTRACT_CONSTRUCT_READS
// already wrote out: this process runs in a python container, with no samtools
// in it.

process READ_ARCHITECTURE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(sam), path(construct_fasta)
    path(script)

    output:
    tuple val(meta), path("${meta.id}.architecture.txt"), emit: report

    script:
    """
    set -euo pipefail

    gzip -cd ${sam} \\
        | python3 ${script} \\
            --construct '${meta.construct_contig}' \\
            --construct-fa ${construct_fasta} \\
        > ${meta.id}.architecture.txt

    head -40 ${meta.id}.architecture.txt
    """

    stub:
    """
    touch ${meta.id}.architecture.txt
    """
}
