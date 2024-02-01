# Dot

Dot is an interactive dot plot viewer for genome-genome alignments.

This is an extended version of Maria Nattestad's original software, which is available
at [github.com/MariaNattestad/dot](https://github.com/MariaNattestad/dot)

Dot is publicly available here: [dot.sandbox.bio](https://dot.sandbox.bio)
And can also be used locally by cloning this repository and simply opening the index.html file in a web browser.

# My modifications

- Create all files required for Dot in one command
- In particular, this script can create Dot-compatible annotation files
- This library provides Python bindings for Mummer4's Nucmer, available here: [github.com/mummer4/mummer](https://github.com/mummer4/mummer)
- Easy installation via pip

# To do:
- [ ] Add optional title string to html
- [ ] Create fna from gbk if possible!

## Installation

This package requires at least `Python 3.9`.

Mummer must be installed.

```bash
pip install git+https://github.com/MrTomRod/dot.git
```

## Usage

**Create all files required for Dot in one command**

```shell
# just the dotplot
dot run --html-out /path/to/dotplot.html --fasta-ref ref.fna --fasta-qry qry.fna
# dotplot + genes
dot run --html-out /path/to/dotplot.html --fasta-ref ref.fna --fasta-qry qry.fna --gbk-ref ref.gbk --gbk-qry qry.gbk
```

This will:

1) run nucmer -> `out.delta`
2) run DotPrep -> `out.coords`, `out.coords.idx`, `out.uniqueAnchorFiltered_l10000.delta.gz`
3) create a reference Dot annotation file -> `out.ref.annotations`
4) create a query Dot annotation file -> `out.qry.annotations`
5) create a standalone html file that can be opened in a browser

To keep the intermediary files, add the argument `--outdir /path/to/outdir`

It is possible to run steps 1 and 2 only, to create a dotplot without genes:

```shell
dot run --outdir /path/to/outdir --fasta-ref ref.fna --fasta-qry qry.fna
```

And to add genes in a second step:

```shell
dot create-annotations --gbk ref.gbk --is-ref True > out.ref.annotations
dot create-annotations --gbk qry.gbk --is-ref False --index out.coords.idx > out.qry.annotations
```

**Dot-compatible annotation files**

The reference annotation file is easy to produce: Simpy parse the fields `ref` (contig ID), `ref_start`, `ref_end`, `name` and `strand` (`+` / `-`)
from the GenBank file.

The query annotation file is a bit trickier: I noticed that nucmer sometimes inverts contigs on the y-axis. As a result, the annotations on this
contig are inverted relative to the dot plot. For this reason, my script reads `out.coords.idx` to learn which contigs are inverted. It then also
inverts the locations of the genes. In addition, my script discardss contigs that are not shown in the dot plot (yes, that can happen).

**Python bindings**

This library gives Python bindings for Nucmer:

```python
from dot import Nucmer

nucmer = Nucmer()

outfile = nucmer.align(
    fasta_ref='ref.fasta', fasta_qry='qry.fasta',
    workdir='/path/to/existing/dir',
    arguments={'--mincluster': 100}
)

print(outfile)  # /path/to/existing/dir/out.delta
```

And also all of the Dot functions described above:

```python
from dot import DotPrep

dotprep = DotPrep()

coords, index = dotprep.run_python(fasta_ref='ref.fasta', fasta_qry='qry.fasta', mincluster=100)
# coords: Dot's out.coords file as string
# index: Dot's out.coords.idx file as string

ref_annos = dotprep.gbk_to_annotation_file(gbk='ref.gbk', is_ref=True)
qry_annos = dotprep.gbk_to_annotation_file(gbk='qry.gbk', is_ref=False, index=index)
# Dot-compatible annotation files as string

ref_annos = dotprep.gbk_to_annotation_dict(gbk='ref.gbk', is_ref=True)
qry_annos = dotprep.gbk_to_annotation_dict(gbk='qry.gbk', is_ref=False, index=index)
# Dot-compatible annotations, as Python dictionaries
```
