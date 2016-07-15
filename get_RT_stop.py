#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO
from argparse import RawTextHelpFormatter
import argparse
import os



#initialize distribution
def init_dist(length_seq):
    distribution = {};
    for transcript in length_seq:
        distribution[transcript] = [];
        for i in range(0, length_seq[transcript]):
            distribution[transcript].append(0);
    return distribution

#Get RTSC counts
def get_rtsc(sf, distribution):
    r = {}
    with open(sf, 'r') as f:
        while(True):
            line = f.readline()
            if not line: break
            tl = line.strip().split('\t');
            t = tl[0].strip()
            s = tl[1].strip()
            if t in distribution:
                distribution[t][int(s)-1] = distribution[t][int(s)-1] + 1;
    r = distribution
    return r
        
#Output RTSC counts           
def write_rtsc(result, distribution):
    with open(result, 'w') as h:
        for t in distribution:
            h.write(t);
            h.write('\n')
            for i in range(0, len(distribution[t])-1):
                h.write(str(distribution[t][i]))
                h.write('\t')
            i = i+1
            h.write(str(distribution[t][i]))
            h.write('\n')
            h.write('\n')

        

def __main__():

    parser = argparse.ArgumentParser(description='Calculate number of RT stops mapped to each nucloetide\n', formatter_class=RawTextHelpFormatter)
    parser.add_argument('mf', metavar="<Mapped_file>",help='The SAM file with only mapped reads and without header')
    parser.add_argument('seq', metavar="<Reference_file>",help='(Transcriptome) reference library used to map the reads (fasta format)')
    parser.add_argument('result', nargs='?', metavar="<RTSC_file>",help='Output RT stop count file (text file) [Default: [Mapped_file name].rtsc]')
    
    args = parser.parse_args()
    
    fasta_file = args.seq
    map_file = args.mf



    syspathrs = os.getcwd()
    map_fn = os.path.basename(map_file)

    if args.result:
        result_file = args.result
    else:
        result_file = map_fn.split(".sam")[0].strip()+".rtsc"

    im = os.path.join(syspathrs, map_fn+"_map_info.txt") #Get the mapping information (intermediate file)

    os.system("cut -f 3,4 "+map_file+" > "+im)
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta');
    length_seq = {}
    for seq in fasta_sequences:
        length_seq[seq.id] = len(str(seq.seq))  #Get length of each transcript

    di = init_dist(length_seq)
    
    du = get_rtsc(im, di)
    
    write_rtsc(result_file, du)
    
    os.system("rm "+im)  #Remove the intermediate file



if __name__ == "__main__" : __main__()



        
        
  






