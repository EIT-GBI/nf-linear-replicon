# Citations

## Pipeline

- **Nextflow** — Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C.
  *Nextflow enables reproducible computational workflows.* Nat Biotechnol. 2017.
  doi: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

## Tools

- **minimap2** — Li H. *Minimap2: pairwise alignment for nucleotide sequences.*
  Bioinformatics. 2018. doi: [10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)
- **samtools / htslib** — Danecek P, Bonfield JK, Liddle J, et al.
  *Twelve years of SAMtools and BCFtools.* GigaScience. 2021.
  doi: [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)
- **dorado** — Oxford Nanopore Technologies. https://github.com/nanoporetech/dorado
  (optional aligner path)

## Origin

The logic in `modules/local/*` is a port of the `wei_linear_replicon` cluster
pipeline (`prep_ref.sh`, `align.sh`, `junction_positions.sh`, `junctions.sh`,
`architecture.py`), including its thresholds and its defences against the
failure modes documented there.
