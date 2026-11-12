// Step 3b: junction recurrence matrix across all samples.
//
// A real integration recurs at one locus across independent libraries of the
// same clone; ONT chimeras scatter and do not. Breakpoints are pooled across
// samples, clustered by position (params.matrix.cluster_tolerance) and counted
// as UNIQUE reads per sample.

process JUNCTION_MATRIX {

    input:
    path(split_bams)
    path(names)

    output:
    path("junction_matrix.csv"),  emit: csv
    path("junction_summary.csv"), emit: summary

    script:
    def tolerance   = params.matrix.cluster_tolerance
    def min_reads   = params.matrix.min_total_reads
    def host_regex  = params.host_contig_prefixes.join('|')
    """
    set -euo pipefail

    : > breakpoints.tsv
    : > summary.tsv
    : > samples.txt

    for bam in \$(ls *.split.bam | sort); do
        id=\$(basename "\$bam" .split.bam)
        echo "\$id" >> samples.txt

        # The construct is the one contig that is not a host accession.
        construct=\$(samtools view -H "\$bam" \\
            | awk '\$1 == "@SQ" { sub("SN:", "", \$2); if (\$2 !~ /^(${host_regex})/) print \$2 }' | head -1)
        [ -n "\$construct" ] || { echo "\$id: could not identify the construct contig" >&2; exit 1; }

        # Reads with both a construct and a genome alignment. -F 0x104 keeps
        # supplementary (0x800): the genome side is usually the supplementary.
        samtools view -F 0x104 "\$bam" \\
            | awk -v construct="\$construct" '
                { if (\$3 == construct) on_construct[\$1] = 1; else on_genome[\$1] = 1 }
                END { for (read in on_genome) if (read in on_construct) print read }' > keep.txt

        samtools view -F 0x104 "\$bam" \\
            | awk -v construct="\$construct" -v id="\$id" '
                NR == FNR { keep[\$1]; next }
                \$3 != construct && (\$1 in keep) { print \$3 "\\t" \$4 "\\t" id "\\t" \$1 "\\t" \$5 }' keep.txt - \\
            >> breakpoints.tsv

        printf "%s\\t%s\\t%s\\t%s\\n" "\$id" "\${id%%__*}" \\
            "\$(wc -l < "\$id.names.txt" | tr -d ' ')" "\$(wc -l < keep.txt | tr -d ' ')" >> summary.tsv
    done

    { echo "sample_id,sample,construct_reads,junction_reads,junction_rate_pct"
      awk -F'\\t' '{ printf "%s,%s,%s,%s,%.1f\\n", \$1, \$2, \$3, \$4, (\$3 ? 100 * \$4 / \$3 : 0) }' summary.tsv
    } > junction_summary.csv

    # Pool breakpoints across samples, cluster by position, count UNIQUE reads.
    sort -k1,1 -k2,2n breakpoints.tsv \\
        | awk -v tolerance=${tolerance} 'BEGIN { OFS = "\\t" }
            { if (\$1 != previous_contig || \$2 - previous_position > tolerance) {
                  cluster++ ; anchor[cluster] = \$1 ":" \$2
              }
              previous_contig = \$1 ; previous_position = \$2
              key = cluster "\\t" \$3 "\\t" \$4
              if (!(key in seen)) { seen[key] = 1 ; reads[cluster "\\t" \$3]++ }
              locus[cluster] = anchor[cluster] }
            END { for (k in reads) { split(k, part, "\\t"); print part[1], locus[part[1]], part[2], reads[k] } }' \\
        > counts.tsv

    samples=\$(paste -sd, samples.txt)
    n_samples=\$(wc -l < samples.txt | tr -d ' ')

    # Pivot: one row per locus, one column per sample.
    { printf "contig,position,%s,total_reads,n_samples\\n" "\$samples"
      awk -F'\\t' -v OFS=, -v columns="\$samples" -v minimum=${min_reads} '
        BEGIN { n = split(columns, sample, ",") }
        { matrix[\$1 "\\t" \$3] = \$4 ; locus[\$1] = \$2 ; total[\$1] += \$4 }
        END {
            for (cluster in locus) {
                if (total[cluster] < minimum) continue
                split(locus[cluster], where, ":")
                line = where[1] "," where[2] ; shared = 0
                for (i = 1; i <= n; i++) {
                    value = matrix[cluster "\\t" sample[i]] + 0
                    line = line "," value
                    if (value) shared++
                }
                print line "," total[cluster] "," shared
            }
        }' counts.tsv \\
      | sort -t, -k\$(( n_samples + 4 )),\$(( n_samples + 4 ))nr \\
             -k\$(( n_samples + 3 )),\$(( n_samples + 3 ))nr
    } > junction_matrix.csv

    cat junction_summary.csv
    echo
    echo "\$(( \$(wc -l < junction_matrix.csv) - 1 )) loci with >= ${min_reads} read(s) in junction_matrix.csv"
    """

    stub:
    """
    echo "contig,position,total_reads,n_samples" > junction_matrix.csv
    echo "sample_id,sample,construct_reads,junction_reads,junction_rate_pct" > junction_summary.csv
    """
}
