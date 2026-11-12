// Step 3a, part 2: the per-sample junction tables, pooled into one report.
//
// START WITH THE DISTRIBUTION TABLE. Accidental joins made during library prep
// land in random places, so they give ONE read per position. A real integration
// is crossed by several reads at the same place. Comparing a control against a
// sample here is how you pick the cutoff -- not a fixed threshold.

process JUNCTION_REPORT {

    input:
    path(per_sample_csv)

    output:
    path("junction_positions.csv"),              emit: csv
    path("junction_positions_distribution.csv"), emit: distribution
    path("junction_positions_candidates.csv"),   emit: candidates

    script:
    def min_reads     = params.junction.min_reads_to_shortlist
    def max_to_print  = params.junction.max_rows_to_print
    """
    set -euo pipefail

    # A sample with no junctions is a result in itself, so it must still appear
    # in the distribution table rather than silently vanishing. The sample list
    # comes from the file names, which exist for every sample either way.
    for f in *.junction_positions.csv; do
        basename "\$f" .junction_positions.csv
    done | sort > samples.txt

    head -1 \$(ls *.junction_positions.csv | head -1) > junction_positions.csv
    for f in \$(ls *.junction_positions.csv | sort); do
        tail -n +2 "\$f"
    done >> junction_positions.csv

    found=\$(( \$(wc -l < junction_positions.csv) - 1 ))

    # --- 1. read support per position ------------------------------------
    {
        echo "sample,positions_with_1_read,with_2_reads,with_3_reads,with_4_reads,with_5_reads,with_6_reads,with_7_or_more"
        tail -n +2 junction_positions.csv \\
            | awk -F, -v sample_list=samples.txt '
                BEGIN { OFS = ","
                        while ((getline line < sample_list) > 0) sample[line] = 1 }
                { count[\$1 SUBSEP \$4]++ ; sample[\$1] = 1 }
                END {
                    for (s in sample) {
                        line = s
                        for (r = 1; r <= 6; r++) line = line OFS (count[s SUBSEP r] + 0)
                        many = 0
                        for (k in count) {
                            split(k, part, SUBSEP)
                            if (part[1] == s && part[2] + 0 >= 7) many += count[k]
                        }
                        print line OFS many
                    }
                }' \\
            | sort -t, -k1,1
    } > junction_positions_distribution.csv

    # --- 2. the shortlist -------------------------------------------------
    {
        head -1 junction_positions.csv
        # by read count, then sequence on each side: the most convincing
        # junctions come top whether or not they have more than one read
        tail -n +2 junction_positions.csv \\
            | awk -F, -v minimum=${min_reads} '\$4 >= minimum' \\
            | sort -t, -k4,4nr -k6,6nr -k9,9nr
    } > junction_positions_candidates.csv

    echo "=== positions with 1 read, 2 reads, and so on ==="
    cat junction_positions_distribution.csv
    echo
    echo "=== junction positions, best supported first (top ${max_to_print}) ==="
    head -n \$(( ${max_to_print} + 1 )) junction_positions_candidates.csv
    echo
    echo "\$found positions in junction_positions.csv"
    """

    stub:
    """
    echo "sample,genome_contig,genome_position,reads,best_mapping_quality,genome_bp,construct_from,construct_to,construct_bp,genome_bp_each,construct_bp_each" > junction_positions.csv
    cp junction_positions.csv junction_positions_candidates.csv
    echo "sample,positions_with_1_read,with_2_reads,with_3_reads,with_4_reads,with_5_reads,with_6_reads,with_7_or_more" > junction_positions_distribution.csv
    """
}
