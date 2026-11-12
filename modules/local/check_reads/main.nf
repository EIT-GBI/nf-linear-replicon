// Fail on a corrupt read file before the alignment starts, not hours into it.
//
// A truncated .fastq.gz desynchronises the SAM stream and samtools aborts much
// later with the misleading "CIGAR and query sequence are of different length".
// For BAMs, `samtools quickcheck` reads each header and verifies the EOF block
// -- a stricter check than gzip -t, and the reason a truncated bam_pass file is
// caught here.
//
// The reads are passed straight through, so this process is a gate in the chain
// rather than a side branch that cannot block anything.

process CHECK_READS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads, includeInputs: true), emit: reads

    script:
    """
    set -euo pipefail

    for f in ${reads}; do
        case "\$f" in
            *.bam)
                # -u: bam_pass is UNALIGNED and has no @SQ lines. Without this
                # every file is reported "had no targets in header".
                samtools quickcheck -v -u "\$f" || { echo "CORRUPT bam: \$f" >&2; exit 1; }
                ;;
            *.gz)
                gzip -t "\$f" || { echo "CORRUPT fastq: \$f" >&2; exit 1; }
                ;;
        esac
    done

    echo "${meta.id}: ${reads instanceof List ? reads.size() : 1} ${meta.read_type} file(s) OK"
    """

    stub:
    """
    true
    """
}
