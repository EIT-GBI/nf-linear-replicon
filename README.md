# nf-linear-replicon

Finds where a linear replicon construct has integrated into a host genome, from
Nanopore reads. It is the `wei_linear_replicon` cluster pipeline
(`prep_ref` → `align` → `analyse`) rewritten as Nextflow, with the same logic,
the same thresholds and the same output files.

**The whole idea in one line:** align each sample against the genome **and** the
construct *together*. A read that crosses an integration junction then gets a
primary alignment on one sequence and a supplementary alignment on the other,
and minimap2 records the pairing in the `SA` tag. Everything below is read off
those pairs.

---

## What the pipeline does

```mermaid
flowchart TD
    A[samplesheet.csv] --> B[BUILD_REFERENCE<br/><i>genome + construct in one FASTA</i>]
    A --> C[PREPARE_READS<br/><i>quickcheck / gzip -t, then one file per sample</i>]
    B --> D[CONSTRUCT_HOMOLOGY<br/><i>construct vs genome: shared sequence</i>]
    B --> E[SAMTOOLS_FAIDX]

    C --> F{alignment.tool}
    F -->|dorado, default| G[DORADO_ALIGNER<br/><i>reads unaligned BAM directly</i>]
    F -->|minimap2| H[MINIMAP2_ALIGN<br/><i>-Y -L --MD --secondary=no</i>]
    B --> G
    B --> H

    G --> I[sorted BAM + index]
    H --> I
    I --> J[SAMTOOLS_FLAGSTAT]
    I --> K[EXTRACT_CONSTRUCT_READS<br/><i>construct_reads.bam + split.bam</i>]
    K --> L[JUNCTION_POSITIONS<br/><i>per sample</i>]
    L --> M[JUNCTION_REPORT<br/><i>distribution + candidates</i>]
    K --> N[JUNCTION_MATRIX<br/><i>recurrence across samples</i>]
    K --> O[READ_ARCHITECTURE<br/><i>optional, per-read structure</i>]
```

| Stage | What it is for |
|---|---|
| `BUILD_REFERENCE` | one combined genome+construct FASTA per pair, and the name of the construct contig |
| `CONSTRUCT_HOMOLOGY` | every alignment between construct and genome. Shared sequence makes fake junctions at a **fixed construct coordinate** in every sample — read this before believing a result |
| `PREPARE_READS` | checks every read file (`samtools quickcheck` / `gzip -t`) — a truncated file otherwise kills samtools hours later with a misleading error — and hands the aligner one file. One file in a folder is just renamed; several are joined with `samtools cat`, a BGZF block copy |
| `DORADO_ALIGNER` / `MINIMAP2_ALIGN` | dorado reads unaligned BAM directly, so `bam_pass` needs no FASTQ round trip. minimap2 adds `-Y`, which keeps the split pieces soft-clipped so the sequence either side of a junction stays in the record |
| `EXTRACT_CONSTRUCT_READS` | the two small BAMs: every alignment of every read touching the construct, and of those, the reads broken across two places |
| `JUNCTION_POSITIONS` | per sample: every genome position where the construct joins the genome, with read support |
| `JUNCTION_REPORT` | those tables pooled: the distribution table and the shortlist |
| `JUNCTION_MATRIX` | the same junctions pooled and clustered across samples — recurrence is the strongest evidence there is |
| `READ_ARCHITECTURE` | optional: per-read segment maps, independent-molecule counts, construct–construct junctions |

---

## Run it

Full step-by-step guide, written for someone who has not used Nextflow:
**[docs/running-the-pipeline.md](docs/running-the-pipeline.md)**. The short
version, on the cluster:

```bash
# 1. a run folder for your settings and results, and a SEPARATE work folder for
#    Nextflow's large temporary data, so you can delete it in one go later.
#    Keep outdir inside the run folder: results are published as hard links
#    into the work folder, and a hard link cannot cross filesystems.
DATASET_NAME=20260928
RUN_DIR=/mnt/lustre/users/$USER/data/wei_linear_replicon/$DATASET_NAME
export NXF_WORK=/mnt/lustre/users/$USER/nf-work/$DATASET_NAME
mkdir -p $RUN_DIR $NXF_WORK
cd $RUN_DIR

# 2. your samplesheet (see below) goes here, then fetch the settings file and
#    edit the three paths at the top of it
curl -O https://raw.githubusercontent.com/EIT-GBI/nf-linear-replicon/main/params.cluster.yaml
nano params.cluster.yaml

# 3. submit. One Slurm job, which is Nextflow itself; it submits the rest.
module load nextflow
sbatch -J nf-linear-replicon -p cpu -t 2-00:00:00 \
  --wrap='nextflow run EIT-GBI/nf-linear-replicon -latest -params-file params.cluster.yaml -profile cluster -resume'
```

`-resume` reuses everything already finished, so changing a threshold re-runs
the analysis in minutes without realigning anything.

On a laptop, the same shape — a folder, the params file edited in it, and a work
folder of your choosing:

```bash
mkdir -p ~/my-run && cd ~/my-run
curl -O https://raw.githubusercontent.com/EIT-GBI/nf-linear-replicon/main/params.cluster.yaml
nano params.cluster.yaml
export NXF_WORK=~/nf-work/my-run
nextflow run EIT-GBI/nf-linear-replicon -latest -params-file params.cluster.yaml -profile docker -resume
```

Or try the bundled example first, which has an integration planted at a known
position:

```bash
nextflow run . -profile test,docker     # or test,local for no containers at all
```

`-profile local` runs without any container, using the tools on your `$PATH`.
`-profile docker` / `-profile apptainer` use the pinned module images.

### The samplesheet

```csv
sample,readset,barcode,reference_genome,additional_reference,fastq_path
LWm010,20260928,barcode01,CHOWT_genome,LP247_maxi,/path/to/run/fastq_pass/barcode01
LWm010,20261023,barcode04,CHOWT_genome,LP247_maxi,/path/to/other_run/bam_pass/barcode04
CONTROL,20260928,barcode09,CHOWT_genome,LP247_maxi,/path/to/run/fastq_pass/barcode09
```

| column | |
|---|---|
| `sample` | the clone. Repeats across readsets on purpose |
| `readset`, `barcode` | optional, but outputs are keyed on `sample__readset__barcode`, so without them two libraries of one clone collide |
| `reference_genome`, `additional_reference` | file names under `--reference_dir`. The extension may be left off: `.fa`, `.fasta`, `.fna` and their `.gz` forms are all tried |
| `fastq_path` | a **folder**, searched recursively. A relative path is taken relative to the samplesheet |

**Headers are matched loosely**, as the cluster scripts matched them:
case-insensitively, trimmed, CRLF- and BOM-safe, and **by prefix** — so
`Sample`, `ReadSet` and `additional_reference_sequence` all work, as do
`additional_reference` and `additional_reference_name`. A column that cannot be
matched is a hard error listing the headers that were found.

**Sample names are reduced to a safe token for file names**: instrument
samplesheets carry things like `CHOWT genome +LP PRD1 maxi (0.1%)`, which would
otherwise reach file names and command lines. That becomes
`CHOWT_genome_LP_PRD1_maxi_0.1`; the original is kept in the sample metadata.

**`fastq_path` may point at `bam_pass` instead of `fastq_pass`.** Dorado writes
unaligned BAM by default now. The column name does not change: FASTQs are looked
for first and `*.bam` is the fallback, and FASTQ wins if a folder holds both.
On the BAM route `samtools fastq -T MM,ML` plus minimap2 `-y` carry the
modified-base tags through, and `samtools quickcheck` replaces the `gzip -t`
precheck.

### Choosing the aligner

`params.alignment.tool` is `dorado` (default) or `minimap2`.

**dorado** is `dorado aligner`, which is minimap2 inside. It reads unaligned BAM
directly, so a `bam_pass` run never becomes FASTQ and the `MM`/`ML`
modified-base tags are never touched. It takes one reads file, which
`PREPARE_READS` has already seen to. It runs on CPU — only basecalling needs a
GPU — so the cluster profile keeps it on the `cpu` queue. It also carries its
own `publishDir`, being a submodule, so on this route the BAMs are published by
the module rather than by this pipeline's `output` block.

**minimap2** is the alternative, and the path the explicit flags are written
for: `-Y -L --MD --secondary=no`. It cannot read BAM, so `PREPARE_READS` decodes
a `bam_pass` row to FASTQ first, carrying the tags with `-T MM,ML`. This is the
route the end-to-end tests exercise, because it runs without the dorado image.

Neither choice changes what the junction steps see: they filter with `-F 0x104`,
which drops secondary alignments whatever the aligner emitted.

## Reading the result

In `<outdir>/results/`:

| file | what it shows |
|---|---|
| `junction_positions_distribution.csv` | how many positions have 1 read, 2 reads, … per sample |
| `junction_positions_candidates.csv` | every position, best-supported first |
| `junction_positions.csv` | the same, unsorted |
| `junction_matrix.csv` | one row per locus, one column per sample |
| `junction_summary.csv` | junction reads as a fraction of construct reads |
| `architecture/<id>.architecture.txt` | optional: per-read structure and the loci it implies |

**Start with the distribution table.** Accidental joins made during library prep
land in random places, so they give one read per position. A real integration is
crossed by several reads at the same place. Comparing a control against a sample
there is how you pick the cutoff — not a fixed threshold.

A real integration site has several reads at the same position, long alignments
on both sides (500 bp each is already enforced), recurs in repeat libraries of
the same clone (`junction_matrix.csv`), and is **absent from your controls**.
Any position that also appears in a no-construct or spike-in control is an
artifact, however good it looks.

Two traps:

- Identical per-read lengths at one position (`6152;6152` in the `*_each`
  columns) usually mean the same molecule was sequenced twice — less
  independent support than the read count suggests. The optional architecture
  report counts distinct shear points for exactly this reason.
- `construct_from` equal to `construct_to` is **not** a 1 bp match. It means
  every read broke at the same construct base, which is a good sign. Use
  `construct_bp` for the length.

Then open `subset/*.split.bam` in IGV at the shortlisted positions. Those files
are small; download them rather than the full BAMs.

Also worth a look before interpreting anything: `ref/*.homology.txt` lists every
alignment between construct and genome. Shared sequence produces fake junctions
at a *fixed construct coordinate* in every sample.

### The rest of the output

```
<outdir>/
├── ref/        combined FASTA + .fai, construct-vs-genome homology (PAF + readable)
├── aligned/    <id>.sorted.bam + .bai
├── qc/
│   ├── flagstat/   alignment metrics
│   └── contigs/    which non-construct contigs the split reads hit
├── subset/     <id>.construct_reads.bam, <id>.split.bam (+ indexes), <id>.names.txt
└── results/    the tables above
```

---

## Tuning

Everything is a parameter; the defaults are the cluster pipeline's.

| param | default | |
|---|---|---|
| `junction.min_genome_bp` / `min_construct_bp` | 500 | the main false-positive filter: reads joined by accident, and the short stretches the construct shares with the genome, both leave one **short** side |
| `junction.group_within` | 100 | Nanopore breakpoints wobble a few bases |
| `junction.min_mapping_quality` | 60 | below this is repetitive sequence, not trustworthy |
| `junction.min_reads_to_shortlist` | 1 | 1 keeps the single-read background you judge against |
| `matrix.cluster_tolerance` | 2000 | breakpoints within this distance are one locus |
| `matrix.min_total_reads` | 2 | drop clusters with fewer reads than this in total |
| `host_contig_prefixes` | `['NC_', 'NW_']` | how the construct contig is told apart from the host's |
| `architecture.enabled` | false | the per-read architecture report |

`host_contig_prefixes` matters more than it looks: everything downstream
identifies the construct as *the contig that is not a host accession*. A
construct whose FASTA header starts with one of these is a hard error in
`BUILD_REFERENCE`, because it would otherwise be taken for a chromosome and
silently yield zero junctions.

---

## Testing

```bash
nf-test test                      # laptop: tools from $PATH
nf-test test --profile test,docker
nextflow run . -profile test,local
```

The test data is in `assets/test-data`, built from a fixed seed by
`scripts/make_test_data.py`, and **its truth is known**: a 75 kb three-contig
genome, a 6 kb construct, one integration planted at `NC_900001.1:20000` joined
at construct base 1200, crossed by 6 reads in one library of `SAMPLE1` and 5 in
the other, absent from `CONTROL`, with one-off chimeras scattered through every
sample as background. `assets/test-data/truth.json` records it.

So the tests assert that the pipeline **found the integration that was
planted** — not that a checksum did not move. They also check that the control
stays clean, that raising `min_construct_bp` above the planted junction removes
it (the filter is really being applied), that unaligned BAM input gives the same
answer as FASTQ, and that the dorado route skips the conversion step.

---

## Layout

```
main.nf                     the whole workflow, top to bottom
nextflow.config             parameters, containers, resources, profiles
modules/samtools|minimap2|dorado   nf-mod-* git submodules (do not edit)
modules/local/              the steps with no module of their own
scripts/make_test_data.py   rebuilds assets/test-data from a fixed seed
tests/main.nf.test          nf-test suite
docs/running-the-pipeline.md   step-by-step guide for running it
```

`modules/local/*` are ports of the cluster shell scripts, one process each, with
their awk kept intact so the two can be compared line by line.

## What is not here

`replicon_reads.py` and `plots.py` from the cluster repo are exploratory tools —
whole-read re-mapping, palindrome and telomere-tract analysis, figures — run by
hand on `subset/` once you have a shortlist. They are not part of the pipeline
and are not bundled; use them from the `wei_linear_replicon` repo against the
published `subset/` folder.
