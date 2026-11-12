#!/usr/bin/env python3
"""Build the tiny test dataset in assets/test-data/, from a fixed seed.

    python3 scripts/make_test_data.py [--outdir assets/test-data]

Nothing here is downloaded and nothing is random between runs: the same seed
gives the same files, byte for byte. What matters about this data is that its
truth is known and written down in assets/test-data/truth.json, so the tests
can assert that the pipeline found the integration that was planted, rather
than that a checksum did not move.

What is planted
---------------
  genome      three contigs named like a RefSeq assembly: NC_900001.1 (40 kb),
              NC_900002.1 (25 kb), NW_900003.1 (10 kb). The NC_/NW_ prefixes
              matter: everything downstream identifies the construct as the one
              contig that is NOT a host accession.
  construct   LP247_mini, 6 kb.
  samples     SAMPLE1 sequenced twice (rep1, rep2) carries ONE integration, at
              NC_900001.1:20000, crossed by several reads in both libraries --
              fixed on the junction side, scattered on the shear side, which is
              what separates a real site from a chimera.
              CONTROL carries the construct but no integration.
  background  every sample gets a few one-off construct/genome chimeras at
              random positions, the library-prep artifacts the distribution
              table in the results is there to let you judge against.

Reads are error-prone (substitutions only) so that mapping is realistic, but
long and unique enough to map at MAPQ 60.
"""
import argparse
import gzip
import json
import random
import shutil
import subprocess
from pathlib import Path

SEED = 20260928

GENOME_CONTIGS = [('NC_900001.1', 40000), ('NC_900002.1', 25000), ('NW_900003.1', 10000)]
CONSTRUCT_NAME = 'LP247_mini'
CONSTRUCT_LEN = 6000

# The one real integration: construct position 1200 joined to genome position
# 20000 on NC_900001.1.
INTEGRATION = {'contig': 'NC_900001.1', 'position': 20000, 'construct_position': 1200}

ERROR_RATE = 0.015
BASES = 'ACGT'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def random_sequence(rng, n):
    return ''.join(rng.choice(BASES) for _ in range(n))


def mutate(rng, seq):
    """Substitutions only: enough to look like ONT, not enough to break mapping."""
    out = list(seq)
    for i, base in enumerate(out):
        if rng.random() < ERROR_RATE:
            out[i] = rng.choice([b for b in BASES if b != base])
    return ''.join(out)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def write_fasta(path, records):
    with open(path, 'w') as fh:
        for name, seq in records:
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + '\n')


def write_fastq(path, reads):
    path.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 so regenerating gives byte-identical files
    with open(path, 'wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', compresslevel=9, mtime=0) as fh:
        for name, seq in reads:
            fh.write((f'@{name}\n{seq}\n+\n' + '?' * len(seq) + '\n').encode())


def genome_read(rng, genome, name):
    contig, seq = rng.choice(list(genome.items()))
    length = rng.randint(2000, 5000)
    start = rng.randint(0, len(seq) - length)
    piece = seq[start:start + length]
    if rng.random() < 0.5:
        piece = revcomp(piece)
    return (name, mutate(rng, piece))


def construct_read(rng, construct, name):
    length = rng.randint(1500, 4000)
    start = rng.randint(0, len(construct) - length)
    piece = construct[start:start + length]
    if rng.random() < 0.5:
        piece = revcomp(piece)
    return (name, mutate(rng, piece))


def junction_read(rng, genome, construct, name, contig, position, construct_position):
    """A read that crosses a construct -> genome junction.

    The junction coordinates are FIXED; how far the read runs on either side is
    not. That is the signature a real integration site leaves, and the shear
    scatter is what lets the architecture report count independent molecules.
    """
    construct_bp = rng.randint(900, construct_position)
    genome_bp = rng.randint(1200, 3000)
    construct_side = construct[construct_position - construct_bp:construct_position]
    genome_side = genome[contig][position:position + genome_bp]
    read = mutate(rng, construct_side + genome_side)
    if rng.random() < 0.5:
        read = revcomp(read)
    return (name, read)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', default='assets/test-data')
    args = parser.parse_args()

    out = Path(args.outdir)
    if out.exists():
        shutil.rmtree(out)
    (out / 'references').mkdir(parents=True)

    rng = random.Random(SEED)

    genome = {name: random_sequence(rng, length) for name, length in GENOME_CONTIGS}
    construct = random_sequence(rng, CONSTRUCT_LEN)

    write_fasta(out / 'references' / 'cho_mini.fa', list(genome.items()))
    write_fasta(out / 'references' / 'LP247_mini.fa', [(CONSTRUCT_NAME, construct)])

    truth = {
        'integration': INTEGRATION,
        'construct_contig': CONSTRUCT_NAME,
        'samples': {},
    }

    # sample id -> how many reads cross the planted integration
    plan = {'SAMPLE1__rep1__barcode01': 6, 'SAMPLE1__rep2__barcode02': 5, 'CONTROL__rep1__barcode03': 0}

    for sample_id, n_junction in plan.items():
        reads = []
        for i in range(40):
            reads.append(genome_read(rng, genome, f'{sample_id}_genome_{i}'))
        for i in range(25):
            reads.append(construct_read(rng, construct, f'{sample_id}_construct_{i}'))
        for i in range(n_junction):
            reads.append(junction_read(rng, genome, construct, f'{sample_id}_junction_{i}',
                                       INTEGRATION['contig'], INTEGRATION['position'],
                                       INTEGRATION['construct_position']))

        # One-off chimeras: the library-prep background. Each lands somewhere
        # different, so each gives a single-read position.
        chimeras = []
        for i in range(3):
            contig = rng.choice(list(genome))
            position = rng.randint(1000, len(genome[contig]) - 4000)
            construct_position = rng.randint(1500, CONSTRUCT_LEN - 500)
            reads.append(junction_read(rng, genome, construct, f'{sample_id}_chimera_{i}',
                                       contig, position, construct_position))
            chimeras.append({'contig': contig, 'position': position})

        rng.shuffle(reads)
        write_fastq(out / 'fastq' / sample_id / 'reads.fastq.gz', reads)
        truth['samples'][sample_id] = {'junction_reads': n_junction, 'chimeras': chimeras}

    # The same reads as an unaligned BAM, so the bam_pass route (what dorado
    # writes by default now) is covered by a test as well.
    if shutil.which('samtools'):
        sample_id = 'SAMPLE1__rep1__barcode01'
        ubam = out / 'ubam' / sample_id
        ubam.mkdir(parents=True)
        subprocess.run(
            ['samtools', 'import', '-0', str(out / 'fastq' / sample_id / 'reads.fastq.gz'),
             '-o', str(ubam / 'reads.unaligned.bam')],
            check=True)
    else:
        print('samtools not on PATH: skipping the unaligned BAM')

    (out / 'truth.json').write_text(json.dumps(truth, indent=2) + '\n')
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
