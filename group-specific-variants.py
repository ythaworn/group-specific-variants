#!/usr/bin/env python3
#
# group-specific-variants.py
#
# Copyright (c) 2021 Yuttapong Thawornwattana
#
# Requirements:
# (1) python packages:
#     - numpy
#     - pyfaidx
#
# (2) inputs
#     - directory to grouplist files, each containing a list of sample names
#       that belong to the same clade, with filename being the clade name


import os
import argparse
import timeit

import numpy as np
from pyfaidx import Fasta


# mapping for nucleotide alignment chatacters
nt2int = {'T': 0, 'C': 1, 'A': 2, 'G': 3, '-': 4, 'N': 5, '?': 6}
int2nt = {0: 'T', 1: 'C', 2: 'A', 3: 'G', 4: '-', 5: 'N', 6: '?'}


# read file and output a list
def scan(filename):
    out = []
    with open(filename, 'r') as f:
        for line in f:
            for l in line.split():
                out.append(l)
    return out


def get_aln(aln_file):
    assert os.path.isfile(aln_file), "error: alignment file doesn't exist"

    aln = Fasta(aln_file)
    sq_name = sorted(aln.keys())


    # convert aln into matrix
    num_sq = len(sq_name)
    num_site = len(aln[sq_name[0]])
    aln_np = np.zeros((num_sq, num_site), dtype=int)

    for i in range(num_sq):
        aln_np[i, :] = [nt2int[x] for x in list(aln[sq_name[i]][:].seq)]

    del aln  # remove aln to free up some memory
    return aln_np, sq_name, num_sq, num_site


def get_ref_sq(ref_file):
    assert os.path.isfile(ref_file), "error: reference file doesn't exist"
    ref = Fasta(ref_file)
    ref_name = ref.keys()

    # sanity check: should be one sq
    assert len(ref_name) == 1, "error: reference file contains more than one sequence"

    ref_sq = ref[list(ref_name)[0]]
    del ref, ref_name
    return ref_sq


def main(args):
    # check if input arguments are valid
    if not os.path.isdir(args.group_dir):
        sys.exit('ERROR: invalid group-list dir: ' + args.group_dir)

    if not os.path.isfile(args.aln):
        sys.exit('ERROR: invalid aln file: ' + args.aln)

    if not os.path.isfile(args.ref):
        sys.exit('ERROR: invalid ref file: ' + args.ref)

    if not os.path.isfile(args.pos):
        sys.exit('ERROR: invalid pos file: ' + args.pos)

    if args.method not in [1,2,3,4]:
        sys.exit('ERROR: invalid method: ' + args.method)

    # create out_dir
    os.makedirs(args.out_dir, exist_ok=True)

    # load aln
    aln_np, sq_name, num_sq, num_site = get_aln(args.aln)
    
    # load group info
    group_list_files = [os.path.join(args.group_dir, f) for f in os.listdir(args.group_dir) if f.endswith('.txt')]
    group_list = {}
    for group_txt in group_list_files:
        group_name = os.path.splitext(os.path.basename(group_txt))[0]  # use string up to first . as group name
        group_list[group_name] = scan(group_txt)
    group_list_keys = sorted(group_list)

    # get ref positions for aln sites
    ref_pos = scan(args.pos)
    ref_pos = [int(x) for x in ref_pos]

    # get ref sq (for output)
    ref_sq = get_ref_sq(args.ref)

    # main: find group-specific variants
    for g in group_list_keys:
        print(str(g))
        count = {}

        # indices of ingroup sqs
        ingr_ind = [sq_name.index(x) for x in group_list[g] if x in sq_name]
        outgr_ind = [x for x in range(num_sq) if x not in ingr_ind]

        # split aln into ingroup and outgroup
        aln_ingr = aln_np[ingr_ind, :]
        aln_outgr = aln_np[outgr_ind, :]
        
        # count base frequencies of columns of aln_ingr
        for j in range(num_site):
            site = aln_ingr[:, j]
            unique = np.unique(site)
            if (len(unique) == 1) and (unique[0] < 4):  # all sqs have the same nt
                count[j] = unique[0]  # record position and type

        # cmp sites in count with outgr
        rm = []
        for site, nt in count.items():
            unique = np.unique(aln_outgr[:, site])
            
            # four criteria
            if ((args.method == 1) and (nt in unique)):  # allow any symbols (incl non-nt) in outgr
                rm.append(site)
            elif ((args.method == 2) and 
                (nt in unique or np.any(np.greater(unique, 3)))):  # don't allow non-nt (symbol > 3)
                rm.append(site)
            elif ((args.method == 3) and 
                (nt in unique or len(unique) > 1)):  # outgr must have the same symbol (can be non-nt)
                rm.append(site)
            elif ((args.method == 4) and 
                (nt in unique or len(unique) > 1 or np.any(np.greater(unique, 3)))):  # outgr must have the same symbol; don't allow gap
                rm.append(site)

        # remove non-group-specific sites
        for i in range(len(rm)):
            del count[rm[i]]
        del rm

        # get ref positions for output
        count_keys = sorted(count)
        pos = [ref_pos[x] for x in count_keys]
        ref_nt = [ref_sq[x - 1].seq for x in pos]  # pos is 1-based index
        gr_nt = [int2nt[count[x]] for x in count_keys]

        # write to file
        outfile = os.path.join(args.out_dir, g + '.txt')
        with open(outfile, 'w') as f:
            hdr = str(len(count)) + '\n'
            f.write(hdr)
            for i in range(len(count)):
                f.write(str(pos[i]) + ' ' + ref_nt[i] + '/' + gr_nt[i] + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract clade-specific variants.')

    # Required positional argument
    parser.add_argument('aln', metavar='aln', type=str,
                        help='directory to vcf files')

    parser.add_argument('ref', metavar='ref', type=str,
                        help='reference genome')

    parser.add_argument('pos', metavar='pos', type=str,
                        help='reference coordinates of alignment sites')

    # Optional arguments
    parser.add_argument('-o', '--out', type=str, 
        dest='out_dir', default='output',
        help='output directory (default: output)')

    parser.add_argument('-g', '--gr_dir', type=str, 
        dest='group_dir', default='group-list',
        help='directory containing clade lists ending with .txt, one line per sample (default: group-list)')

    parser.add_argument('-m', '--method', type=int, 
        dest='method', default=1,
        help='criterion for ourgroup when counting clade-specific SNPs (default: 1)')

    args = parser.parse_args()

    start = timeit.default_timer()
    main(args)
    elapsed = timeit.default_timer() - start
    print("elapsed time: %.2f s" % elapsed)
