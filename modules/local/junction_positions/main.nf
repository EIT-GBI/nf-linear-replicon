// Step 3a: every genome position where the construct joins the genome, with the
// number of supporting reads, for one sample.
//
// Thresholds come from params.junction (see nextflow.config):
//   min_genome_bp / min_construct_bp   the main false-positive filter. Reads
//       joined by accident in library prep, and the short stretches the
//       construct shares with the genome, both leave one SHORT side.
//   group_within      Nanopore breakpoints wobble a few bases
//   min_mapping_quality   below this is repetitive sequence, not trustworthy

process JUNCTION_POSITIONS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(split_bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.junction_positions.csv"), emit: csv

    script:
    def construct       = meta.construct_contig
    def min_genome      = params.junction.min_genome_bp
    def min_construct   = params.junction.min_construct_bp
    def group_within    = params.junction.group_within
    def min_quality     = params.junction.min_mapping_quality
    """
    set -euo pipefail

    echo "sample,genome_contig,genome_position,reads,best_mapping_quality,genome_bp,construct_from,construct_to,construct_bp,genome_bp_each,construct_bp_each" \\
        > ${meta.id}.junction_positions.csv

    # One line per qualifying junction read:
    #   contig position read quality genome_bp construct_start construct_bp
    # -F 0x104 keeps supplementary alignments: on a junction read the construct
    # is usually the primary and the genome side the supplementary, so dropping
    # them would hide half of every read.
    samtools view -F 0x104 ${split_bam} \\
        | awk -v construct='${construct}' \\
              -v min_quality=${min_quality} \\
              -v min_genome=${min_genome} \\
              -v min_construct=${min_construct} '
            # reference bases covered, from the CIGAR: M/D/=/X consume reference
            function aligned_length(cigar,   total, operation, number) {
                total = 0
                while (match(cigar, /^[0-9]+[MIDNSHP=X]/)) {
                    operation = substr(cigar, RSTART + RLENGTH - 1, 1)
                    number    = substr(cigar, RSTART, RLENGTH - 1) + 0
                    if (operation ~ /[MD=X]/) total += number
                    cigar = substr(cigar, RSTART + RLENGTH)
                }
                return total
            }
            \$3 == construct { construct_start[\$1] = \$4; construct_bp[\$1] = aligned_length(\$6); next }

            # Genome side stored, not printed: a read is only a junction once
            # BOTH sides are seen, and the construct alignment may come later.
            \$5 >= min_quality { genome_line[++count] = \$3 " " \$4 " " \$1 " " \$5 " " aligned_length(\$6) }

            END {
                for (i = 1; i <= count; i++) {
                    split(genome_line[i], field, " ")
                    read = field[3]
                    if (!(read in construct_bp))            continue   # never reaches the construct
                    if (field[5]            < min_genome)   continue   # too short to believe
                    if (construct_bp[read]  < min_construct) continue
                    print field[1], field[2], read, field[4], field[5], \\
                          construct_start[read], construct_bp[read]
                }
            }' \\
        | sort -k1,1 -k2,2n > per_read.txt

    # Group nearby positions, one summary line each.
    awk -v sample='${meta.id}' -v group_within=${group_within} '
        function report_group() {
            if (read_count == 0) return
            printf "%s,%s,%d,%d,%d,%d,%d,%d,%d,%s,%s\\n",
                   sample, contig, position, read_count, best_quality,
                   genome_bp_total / read_count, construct_from, construct_to,
                   construct_bp_total / read_count, genome_bp_list, construct_bp_list
        }
        {
            if (\$1 != contig || \$2 - previous > group_within) {
                report_group()
                contig = \$1 ; position = \$2
                read_count = 0 ; best_quality = 0
                genome_bp_total = 0 ; construct_bp_total = 0
                construct_from = "" ; construct_to = ""
                genome_bp_list = "" ; construct_bp_list = ""
                delete already_counted
            }
            previous = \$2

            if (!(\$3 in already_counted)) {          # a read can align here twice
                already_counted[\$3] = 1
                read_count++
                genome_bp_total    += \$5
                construct_bp_total += \$7
                genome_bp_list    = genome_bp_list    (read_count > 1 ? ";" : "") \$5
                construct_bp_list = construct_bp_list (read_count > 1 ? ";" : "") \$7
            }
            if (\$4 > best_quality) best_quality = \$4

            construct_end = \$6 + \$7 - 1              # widest construct stretch here
            if (construct_from == "" || \$6 < construct_from) construct_from = \$6
            if (construct_to   == "" || construct_end > construct_to) construct_to = construct_end
        }
        END { report_group() }
    ' per_read.txt >> ${meta.id}.junction_positions.csv

    echo "${meta.id}: \$(( \$(wc -l < ${meta.id}.junction_positions.csv) - 1 )) junction position(s)"
    """

    stub:
    """
    echo "sample,genome_contig,genome_position,reads,best_mapping_quality,genome_bp,construct_from,construct_to,construct_bp,genome_bp_each,construct_bp_each" \\
        > ${meta.id}.junction_positions.csv
    """
}
