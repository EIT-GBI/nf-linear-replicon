// Sequence the construct shares with the genome produces fake junctions at a
// FIXED construct coordinate in every sample. This is the report to read before
// believing any junction.
//
// Every hit is printed, not just the long ones: an empty PAF and a filtered one
// look identical otherwise.

process CONSTRUCT_HOMOLOGY {
    tag "${tag}"

    input:
    tuple val(tag), path(genome), path(construct)

    output:
    tuple val(tag), path("${tag}.homology.paf"), path("${tag}.homology.txt"), emit: report

    script:
    def args = task.ext.args ?: '-c -x asm20'
    """
    set -euo pipefail

    minimap2 ${args} -t ${task.cpus} ${genome} ${construct} > ${tag}.homology.paf

    {
        hits=\$(wc -l < ${tag}.homology.paf | tr -d ' ')
        echo "--- construct vs genome: \$hits alignment(s) ---"
        if [ "\$hits" -gt 0 ]; then
            echo "construct_start construct_end aln_len matches -> genome locus"
            sort -k11,11nr ${tag}.homology.paf \\
                | awk '{print \$3, \$4, \$11, \$10, "->", \$6":"\$8"-"\$9}'
        else
            echo "none -- no construct/genome sequence sharing to mask downstream"
        fi
    } > ${tag}.homology.txt

    cat ${tag}.homology.txt
    """

    stub:
    """
    touch ${tag}.homology.paf ${tag}.homology.txt
    """
}
