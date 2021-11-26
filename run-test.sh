#!/bin/bash

wd=test
ref=$wd/ref.fasta
aln=$wd/aln.fasta
pos=$wd/pos.txt
group_dir=$wd/group-list
out_dir=$wd/output

python3 group-specific-variants.py $aln $ref $pos -g $group_dir -o $out_dir
