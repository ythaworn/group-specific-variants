# Identify group-specific variants

Extract nucleotide variants specific to a group of samples, usually a clade in a phylogeny, from a sequence alignment.


## Input

1. sequence alignment file (fasta format)
2. reference genome (fasta format)
3. positions of alignment sites in the reference genome (`pos.txt`)
4. directory containing lists of sequence names considered to be in the same group (`group-list/<group>.txt`)


## Usage

`python group-specific-variants.py <aln> <ref> <pos> -g <group_dir> -o <out_dir>`

Inputs:
- `<aln>` is the sequence alignment (fasta format)
- `<ref>` is the reference genome (fasta format)
- `<pos>` is a text file containing the positions of alignment sites in the reference genome, one per line
- `<group_dir>` is the directory containing lists of sequence names considered to be in the same group
- `<out_dir>` is the output directory


## Outputs


Each file in the output directory looks like:

```
56
27199 G/A
75233 C/A
235681 G/A
..
```

First line: number of group-specific SNPs found
Each of the subsequent lines: position in the reference genome and ref/non-ref variants


## Test run

An example is provided in the `test` directory. Uncompress aln.fasta.gz and ref.fasta.gz, and then run `run-test.sh`.


## Requirements

Python3 with numpy and pyfaidx.
