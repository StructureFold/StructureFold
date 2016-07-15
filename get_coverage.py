#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from header.parse_dis_react import *
from Bio import SeqIO
from argparse import RawTextHelpFormatter
import argparse





def __main__():
    parser = argparse.ArgumentParser(description='Calculate coverage on each transcript\n', formatter_class=RawTextHelpFormatter)
    parser.add_argument('rtsc', metavar="<RTSC>",help='input RT stop count file')
    parser.add_argument('seq', metavar="<Reference_file>",help='(Transcriptome) reference library used to map the reads (fasta format)')
    parser.add_argument('result', nargs='?', metavar="<coverage_file>",help='Output coverage file (text file) [Default: [RTSC name]_coverage.txt]')
    
    args = parser.parse_args()

    
    dist_file = args.rtsc
    seq_file = args.seq
    output_file = args.result


    n_m = dist_file.split(".rtsc")[0].strip()
    prefix = n_m+"_"

    if output_file == None:
        output_file = prefix+"coverage.txt"

    
    
    seqs = SeqIO.parse(seq_file, 'fasta')
    dist = parse_dist(dist_file)
    dist = dist[1]

    seq_t = {}
    for seq in seqs:
        seq_t[seq.id] = str(seq.seq)

    with open(output_file, 'w') as h:
        for t in dist:
            l = 0
            s = 0
            for i in range(1, len(dist[t])):
                if seq_t[t][i-1] == 'A' or seq_t[t][i-1] == 'C': #only calculate RT count on A/C for DMS experiments
                    s = s+int(dist[t][i])
                    l = l+1
            h.write(t+'\t'+str(float(s)/l)+'\n') #output two columns (1. transcript id; 2. coverage for the transcript)
            

if __name__ == "__main__" : __main__()








        





