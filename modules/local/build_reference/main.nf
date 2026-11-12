// Step 1: one combined genome+construct reference per pair.
//
// Aligning both together is what makes integration detectable: a junction read
// gets a primary alignment on one sequence and a supplementary on the other,
// and minimap2 records the pairing in the SA tag.

process BUILD_REFERENCE {
    tag "${tag}"

    input:
    tuple val(tag), path(genome), path(construct)

    output:
    tuple val(tag), path("${tag}.fa"),           emit: fasta
    tuple val(tag), path("${tag}.contig"),       emit: contig
    tuple val(tag), path("${tag}.construct.fa"), emit: construct_fasta

    script:
    """
    set -euo pipefail

    # NOT cat: coreutils >=9.0 uses copy_file_range() when stdout is a regular
    # file, and that returns ENODATA on Lustre -- reported against the *input*
    # filename, which looks exactly like data corruption.
    uncompress() {
        if [ "\${1##*.}" = "gz" ]; then zcat "\$1"; else dd if="\$1" bs=4M status=none; fi
    }

    uncompress ${construct} > ${tag}.construct.fa
    { uncompress ${genome}; cat ${tag}.construct.fa; } > ${tag}.fa

    grep -m1 '^>' ${tag}.construct.fa | sed 's/^>//' | awk '{print \$1}' > ${tag}.contig
    construct_contig=\$(cat ${tag}.contig)
    echo "construct contig: \$construct_contig"

    # Everything downstream identifies the construct as the one contig that is
    # not a host accession. A construct named this way would be taken for a
    # chromosome and silently yield zero junctions.
    case "\$construct_contig" in
        ${params.host_contig_prefixes.collect { prefix -> "${prefix}*" }.join('|')})
            echo "ERROR: construct contig '\$construct_contig' looks like a host accession; rename it in ${construct}" >&2
            exit 1
            ;;
    esac

    echo "construct length: \$(awk '/^>/ {next} {n += length(\$0)} END {print n}' ${tag}.construct.fa) bp"
    """

    stub:
    """
    touch ${tag}.fa ${tag}.construct.fa
    echo CONSTRUCT > ${tag}.contig
    """
}
