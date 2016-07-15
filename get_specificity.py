#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from header.parse_dis_react import *
from Bio import SeqIO
from argparse import RawTextHelpFormatter
import argparse


def cal_percent(dic):
    r = {}
    sum_n = 0
    for n in dic:
        sum_n = sum_n+dic[n]

    for n in dic:
        r[n] = float(dic[n])/sum_n
    return r
        



def __main__():
    parser = argparse.ArgumentParser(description='Calculate coverage on each transcript\n', formatter_class=RawTextHelpFormatter)
    parser.add_argument('rtsc', metavar="<RTSC>",help='input RT stop count file')
    parser.add_argument('seq', metavar="<Reference_file>",help='(Transcriptome) reference library used to map the reads (fasta format)')
    parser.add_argument('result', nargs='?', metavar="<specificity_file>",help='Output specificity (text file) [Default: [RTSC name]_specificity.txt]')
    
    args = parser.parse_args()

    
    dist_file = args.rtsc
    seq_file = args.seq
    output_file = args.result


    n_m = dist_file.split(".rtsc")[0].strip()
    prefix = n_m+"_"

    if output_file == None:
        output_file = prefix+"specificity.txt"

    
    
    seqs = SeqIO.parse(seq_file, 'fasta')
    dist = parse_dist(dist_file)
    dist = dist[1]

    seq_t = {}
    for seq in seqs:
        seq_t[seq.id] = str(seq.seq)

    nt = {}

    
    for t in dist:
        l = 0
        s = 0
        for i in range(1, len(dist[t])):
            n = seq_t[t][i-1]
            if n not in nt:
                nt[n] = 0
            nt[n] = nt[n]+int(dist[t][i])

    #print(nt)

    pt = cal_percent(nt)
    

    with open(output_file, 'w') as h:
        for n in pt:
            h.write(n+"\t"+str(pt[n])+"\n") #output two columns (1. Nucleotide; 2. Percentage)
        
            

if __name__ == "__main__" : __main__()








        





