# Running the nf-linear-replicon pipeline

**[nf-linear-replicon](https://github.com/EIT-GBI/nf-linear-replicon)** takes
Nanopore reads from clones that were given a linear replicon construct, and
tells you **where in the host genome the construct integrated** — with the read
support behind each position, and the recurrence of each position across your
libraries.

You will normally run it on the **GBI cluster** (Part 1). You can also run it on
your own computer (Part 2), which is useful for trying it out and for tuning
thresholds on a small dataset. You do not need to know Nextflow.

---

## Contents

- **[Part 1 — Running on the GBI cluster](#part-1--running-on-the-gbi-cluster)**
  - [Before you start](#before-you-start)
  - [Step 1 — Log in](#step-1--log-in)
  - [Step 2 — Make your run folder and your work folder](#step-2--make-your-run-folder-and-your-work-folder)
  - [Step 3 — Put your references in place](#step-3--put-your-references-in-place)
  - [Step 4 — Write your samplesheet](#step-4--write-your-samplesheet)
  - [Step 5 — Download the settings file and edit it](#step-5--download-the-settings-file-and-edit-it)
  - [Step 6 — Start the run](#step-6--start-the-run)
  - [Step 7 — Check how it is going](#step-7--check-how-it-is-going)
  - [Step 8 — Read the results](#step-8--read-the-results)
  - [Step 9 — Clean up and move the results off Lustre](#step-9--clean-up-and-move-the-results-off-lustre)
  - [Running it again](#running-it-again)
- **[Part 2 — Running on your own computer](#part-2--running-on-your-own-computer)**

---

## Words you will see

| Word | What it means here |
| --- | --- |
| **Nextflow** | The program that runs the pipeline: it reads the recipe and runs each step in the right order. |
| **container** | A sealed box holding a tool and everything it needs. The pipeline fetches these, so you never install minimap2 or samtools yourself. |
| **profile** | Which setup to use, chosen with `-profile`. `cluster` means "run on the cluster"; `docker` means "run here with Docker". |
| **params file** | A small text file holding your settings — where your data is, where results go. **This is the file you edit.** |
| **run folder** | Where your settings file and your results live. |
| **work folder** | Where Nextflow keeps its large temporary data. Kept apart so you can delete it in one go. |
| **construct** | Your linear replicon. In the samplesheet it is the `additional_reference` column. |

---
---

# Part 1 — Running on the GBI cluster

## Before you start

You need:

1. **A Sandpit account you can log into** — see the
   [access](https://github.com/EIT-GBI/scientific-computing-docs/blob/main/sections/access.md)
   and [SSH](https://github.com/EIT-GBI/scientific-computing-docs/blob/main/sections/ssh.md)
   guides.
2. **Your reads already on the cluster**, as `fastq_pass` or `bam_pass` folders
   (one folder per barcode). Both work.
3. **Two FASTA files**: the host genome, and your construct.
4. **At least one control sample.** A no-construct or spike-in control is how
   you tell a real integration from library-prep background. The pipeline will
   run without one; interpreting it is guesswork.

Check Nextflow is available:

```bash
module avail nextflow
```

> When you start a run you submit **one** Slurm job: Nextflow itself. It then
> submits more Slurm jobs, one per pipeline step, and waits for them. So
> `squeue` shows one long-running job plus a changing set of short ones. That is
> normal.

---

## Step 1 — Log in

```bash
sft login --team eitoxford
ssh sandpit-tokyo-login
```

---

## Step 2 — Make your run folder and your work folder

Each dataset gets two folders on Lustre, not in your home folder:

- a **run folder**, where your settings file and results go;
- a **work folder**, where Nextflow keeps its large temporary data — the BAMs it
  makes along the way. Keeping it apart means you can delete it in one go in
  [Step 9](#step-9--clean-up-and-move-the-results-off-lustre).

Give your dataset a name — the sequencing run date is a good one:

```bash
DATASET_NAME=20260928
```

Then make both folders and go into the run folder:

```bash
RUN_DIR=/mnt/lustre/users/$USER/data/wei_linear_replicon/$DATASET_NAME
export NXF_WORK=/mnt/lustre/users/$USER/nf-work/$DATASET_NAME
mkdir -p $RUN_DIR $NXF_WORK
cd $RUN_DIR
```

`NXF_WORK` is the setting Nextflow reads to know where its work folder is.
**Set it before you submit, every time**, or Nextflow drops a `work/` folder
inside your run folder instead.

> These shell settings only last until you log out. If you log out before
> [Step 6](#step-6--start-the-run), run both blocks again with the same dataset
> name before you submit.

> **Keep `outdir` inside the run folder** (the settings file has it as
> `./results`). The pipeline publishes results as **hard links** into the work
> folder, so they cost no extra disk space — but a hard link cannot cross
> filesystems, so an `outdir` on a different disk from `NXF_WORK` fails.

**Do not copy reads that are already on Lustre.** Point the samplesheet at
them where they are, or make a folder of shortcuts — a shortcut takes no space:

```bash
mkdir -p reads
ln -s /mnt/lustre/path/to/run/fastq_pass/barcode01 reads/
```

---

## Step 3 — Put your references in place

The pipeline needs one folder holding the host genome **and** each construct,
as separate FASTA files:

```text
/mnt/lustre/users/YOU/references/wei_linear_replicon/
├── CHOWT_genome.fa
├── LP247_maxi.fa
└── Integration_control_W55.fa
```

```bash
REF_DIR=/mnt/lustre/users/$USER/references/wei_linear_replicon
mkdir -p $REF_DIR
```

Again: if a genome is already on Lustre, link it rather than copying it:

```bash
ln -s /mnt/lustre/path/to/CHOWT_genome.fa $REF_DIR/
```

No indexing step to run first — the pipeline builds the combined
genome+construct reference and its index itself, once per genome/construct pair,
and reuses it for every sample that asks for the same pair.

> **Your construct's FASTA header must not start with `NC_` or `NW_`.**
> Everything downstream identifies the construct as *the contig that is not a
> host accession*. The pipeline stops with a clear error if it does, rather than
> treating your construct as a chromosome and reporting zero junctions.

---

## Step 4 — Write your samplesheet

In the run folder, make `samplesheet.csv`:

```csv
sample,readset,barcode,reference_genome,additional_reference,fastq_path
LWm010,20260928,barcode01,CHOWT_genome,LP247_maxi,/mnt/lustre/.../20260928/fastq_pass/barcode01
LWm010,20261023,barcode04,CHOWT_genome,LP247_maxi,/mnt/lustre/.../20261023/bam_pass/barcode04
LWm011,20260928,barcode02,CHOWT_genome,LP247_maxi,/mnt/lustre/.../20260928/fastq_pass/barcode02
CONTROL,20260928,barcode09,CHOWT_genome,LP247_maxi,/mnt/lustre/.../20260928/fastq_pass/barcode09
```

| column | |
| --- | --- |
| `sample` | the clone. It is meant to repeat across readsets — that is what the recurrence matrix is built on |
| `readset`, `barcode` | which sequencing run and barcode. Outputs are named `sample__readset__barcode`, so without these two libraries of one clone overwrite each other |
| `reference_genome`, `additional_reference` | file names in your reference folder. The extension is optional: `.fa`, `.fasta`, `.fna` and their `.gz` forms are all tried |
| `fastq_path` | the **folder** of reads, searched recursively. `fastq_pass` or `bam_pass` — FASTQs are looked for first and `*.bam` is the fallback |

A relative `fastq_path` is taken relative to the samplesheet, so a samplesheet
and its reads can move together.

---

## Step 5 — Download the settings file and edit it

Fetch the ready-made settings file into your run folder and open it:

```bash
curl -O https://raw.githubusercontent.com/EIT-GBI/nf-linear-replicon/main/params.cluster.yaml
nano params.cluster.yaml
```

Change the three paths at the top, replacing `YOU` with your username
(`echo $USER` if unsure):

```yaml
samplesheet:   "/mnt/lustre/users/YOU/data/wei_linear_replicon/20260928/samplesheet.csv"
reference_dir: "/mnt/lustre/users/YOU/references/wei_linear_replicon"
outdir:        "./results"
```

Everything else in the file is already at the defaults the cluster pipeline
used. They are written out rather than left implicit so that each run records
the thresholds it was given. The ones worth knowing:

| Setting | What it does |
| --- | --- |
| `junction: min_genome_bp` / `min_construct_bp` | 500 each. The main false-positive filter: a junction must have a long alignment on **both** sides |
| `junction: min_mapping_quality` | 60. Below this is repetitive sequence |
| `junction: group_within` | 100 bp. Nanopore breakpoints wobble a few bases |
| `matrix: cluster_tolerance` | 2000 bp. Breakpoints this close are called one locus |
| `matrix: min_total_reads` | 2. Clusters with fewer reads in total are dropped from the matrix |
| `architecture: enabled` | `false`. Set `true` for the per-read structure report once you have a shortlist |
| `alignment: tool` | `minimap2`. Leave it unless a sample is a single BAM file |

Save with `Ctrl-O`, `Enter`, then `Ctrl-X`.

---

## Step 6 — Start the run

From your run folder, check the work folder is still set — this should print the
`nf-work` path from Step 2, not an empty line:

```bash
echo $NXF_WORK
```

Then submit:

```bash
module load nextflow
sbatch -J nf-linear-replicon -p cpu -t 2-00:00:00 \
  --wrap='nextflow run EIT-GBI/nf-linear-replicon -latest -params-file params.cluster.yaml -profile cluster -resume'
```

**You should see:**

```text
Submitted batch job 519583
```

Write that number down. The run continues on its own, even if you close your
laptop.

| Part | What it does |
| --- | --- |
| `module load nextflow` | Makes Nextflow available. The job inherits it, and `NXF_WORK`, from your session |
| `-J nf-linear-replicon` | Names the job, so you can spot it in `squeue` |
| `-p cpu` | This job is only Nextflow itself; the pipeline steps ask for their own nodes |
| `-t 2-00:00:00` | Two days. Nextflow lives as long as the whole pipeline |
| `-latest` | Fetch the newest version of the pipeline. Swap for `-r v1.0.0` to pin a version |
| `-profile cluster` | Use Slurm and the cluster's containers |
| `-resume` | Reuse anything already finished. Always safe to include |

---

## Step 7 — Check how it is going

```bash
squeue -u $USER
tail -30 slurm-519583.out        # your own job number
```

A step showing `PD` is queued, not broken. Expect the alignment steps to be the
long ones — the CHO genome index is built once per genome/construct pair, then
reused.

> Please check this way, rather than leaving a command running that follows the
> log: the login node is shared, and a run continues perfectly well without
> anyone watching it.

Once it has left `squeue`:

```bash
sacct -j 519583 --format=JobID,JobName,State,Elapsed
```

To stop a run, cancel the Nextflow job — it cleans up the rest:

```bash
scancel 519583
```

---

## Step 8 — Read the results

```text
results/
├── ref/        the combined reference, and the construct-vs-genome homology report
├── aligned/    one sorted BAM per library
├── qc/         flagstat, and which contigs the split reads hit
├── subset/     the small BAMs: construct reads, and the split ones
└── results/    the tables
```

**Start with `results/results/junction_positions_distribution.csv`.** It counts
how many positions have 1 read, 2 reads, and so on, per sample. Accidental joins
made during library prep land in random places, so they give one read per
position; a real integration is crossed by several reads at the same place.
Comparing a control against a sample here is how you pick your cutoff — there is
no fixed threshold that is right for every run.

Then `junction_positions_candidates.csv` (every position, best supported first)
and `junction_matrix.csv` (one row per locus, one column per library). A real
integration site:

- has several reads at the same position,
- has long alignments on both sides (500 bp each is already enforced),
- recurs in repeat libraries of the same clone,
- and is **absent from your controls**. Any position that also appears in a
  control is an artifact, however good it looks.

Two traps:

- Identical per-read lengths at one position (`6152;6152` in the `*_each`
  columns) usually mean the same molecule was sequenced twice — less independent
  support than the read count suggests.
- `construct_from` equal to `construct_to` is **not** a 1 bp match. It means
  every read broke at the same construct base, which is a good sign. Use
  `construct_bp` for the length.

Before believing anything, read `results/ref/*.homology.txt`: sequence the
construct shares with the genome makes fake junctions at a *fixed construct
coordinate* in every sample.

Finally, download `results/subset/*.split.bam` (small) and open them in IGV at
the positions you shortlisted — not the full BAMs.

---

## Step 9 — Clean up and move the results off Lustre

If you have logged out since Step 2, set the names again first:

```bash
DATASET_NAME=20260928
cd /mnt/lustre/users/$USER/data/wei_linear_replicon/$DATASET_NAME
```

The published results are hard links into the work folder — the same files under
two names, costing space only once — so deleting the work folder leaves them
intact:

```bash
rm -rf /mnt/lustre/users/$USER/nf-work/$DATASET_NAME
```

> Only once you are happy with the run: a future `-resume` then has to start
> from scratch.

Then move the results off Lustre. Lustre is fast but expensive and shared, and
only works as scratch space if finished data moves off promptly:

```bash
module load gbi
gbi data move --detach \
  /mnt/lustre/users/$USER/data/wei_linear_replicon/$DATASET_NAME/results \
  /mnt/user-data/$USER/wei_linear_replicon/$DATASET_NAME/results
```

---

## Running it again

Changing a threshold and re-running is cheap: with `-resume`, the alignments are
reused and only the analysis steps run again, in minutes.

```bash
sbatch -J nf-linear-replicon -p cpu -t 2-00:00:00 \
  --wrap='nextflow run EIT-GBI/nf-linear-replicon -latest -params-file params.cluster.yaml -profile cluster -resume --junction.min_construct_bp 1000'
```

A setting given on the command line like that wins over the params file, which
is handy for a one-off. For anything you want to keep, edit the params file so
the run records it.

Adding samples later: add the rows to `samplesheet.csv` and re-run with
`-resume`. The samples already done are not realigned; the cross-sample tables
are rebuilt with everything.

---
---

# Part 2 — Running on your own computer

Useful for trying the pipeline out, and for tuning thresholds on one small
dataset. Not for a CHO-sized genome: the minimap2 index alone needs ~8 GB.

You need **Docker Desktop** running, plus Java and Nextflow:

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

Check it works on the bundled example, which carries an integration planted at a
known position:

```bash
git clone --recurse-submodules https://github.com/EIT-GBI/nf-linear-replicon
cd nf-linear-replicon
nextflow run . -profile test,docker
cat results-test/results/junction_positions_candidates.csv
```

The top two rows should be `NC_900001.1` at position 20001, with 6 and 5 reads,
in the two `SAMPLE1` libraries — and the control should only ever show
single-read background.

For your own data, the steps are the same as Part 1: make a folder, put the
samplesheet in it, download the params file and edit it, then run with your own
work folder:

```bash
mkdir -p ~/my-run && cd ~/my-run
curl -O https://raw.githubusercontent.com/EIT-GBI/nf-linear-replicon/main/params.cluster.yaml
nano params.cluster.yaml          # the three paths at the top

export NXF_WORK=~/nf-work/my-run
nextflow run EIT-GBI/nf-linear-replicon -latest \
    -params-file params.cluster.yaml -profile docker -resume
```

`-profile local` instead of `docker` uses `minimap2`, `samtools` and `python3`
from your `$PATH`, if you would rather not use containers.
