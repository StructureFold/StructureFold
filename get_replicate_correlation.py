#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from header.parse_dis_react import *
from header.write_file import *
from argparse import RawTextHelpFormatter
import argparse
from Bio import SeqIO


def get_name(f_list):
    r = ""
    for i in range(len(f_list)-1):
        r = r+f_list[i].split(".rtsc")[0].strip()+"_"
    i = i+1
    r = r+f_list[i].split(".rtsc")[0]
    return r

def find_common(d_list):
    r = set()
    for t in d_list[0]:
        flag = 1
        for i in range(1, len(d_list)):
            if t not in d_list[i]:
                flag = 0
        if flag == 1:
            r.add(t)
    return r


def __main__():
    parser = argparse.ArgumentParser(description='Calculate correlation between biological replicate libraries\n', formatter_class=RawTextHelpFormatter)
    parser.add_argument('replicate',help='input RT stop count files (rtsc files) of replicates', nargs="+")
    #parser.add_argument('rep2', metavar="<RTSC_file2>",help='input RT stop count file of replicate 2')
    parser.add_argument('-o', dest = 'result',help='Output correlation file (csv file) [Default: [RTSC name]_correlation.csv]')
    parser.add_argument("-d", "--id_list", dest="id_file", default = None, help='The transcript ID list for calculating RT stops')
    parser.add_argument("-s", "--specificity", dest="spe", default = None, help='Type of nucleotide for calculating RT stops [Default: all type of nucleotides]')
    parser.add_argument("-re", "--reference", dest="ref", default = None, help='Mapping reference file (cDNA library) [Must specify if -s is used]')
    
    args = parser.parse_args()
    
    file_list = args.replicate
    output_file = args.result
    ids = args.id_file
    spe = args.spe
    ref = args.ref

    #print(spe)

    if len(file_list)<2:
        sys.exit("Need at least 2 replicates for correlation!")



    if output_file == None:
        prefix = get_name(file_list)
        output_file = prefix+"_correlation.csv"

    spes = set()
    if spe:
        for i in range(len(spe)):
            spes.add(spe[i])

    idss = set()
    if ids:
        idst = read_t_file(idss)
        for i in range(len(idst)):
            idss.add(idst[i][0])

    if ref:       
        refs = SeqIO.parse(ref, 'fasta')
        cdnas = {}
        for seq in refs:
            cdnas[seq.id] = str(seq.seq)
        

    dist_list = []
    for i in range(len(file_list)):
        dist = parse_dist(file_list[i])
        dist = dist[1]
        dist_list.append(dist)



    r = []
    temp = ['ID', 'Position']
    for i in range(len(file_list)):
        temp.append(file_list[i])
    r.append(temp)
        

    
    t_set = find_common(dist_list)
    

    for t in t_set:
        flag_ids = 0
        if ids:
            if t in idss:
                flag_ids = 1
        else:
            flag_ids = 1
        if flag_ids == 1:
            for i in range(1, len(dist_list[0][t])):
                flag_spe = 0
                if spe:
                    if ref:
                        if t in cdna:
                            n = cdna[t][i-1]
                            if n in spes:
                                flag_spe = 1
                        else:
                            sys.stderr.write(t+' not in reference\n')
                    else:
                        sys.stderr.write('Please specify mapping reference file\n')
                else:
                    flag_spe = 1
                if flag_spe == 1:
                    temp = [t, str(i), dist_list[0][t][i]]
                    for j in range(1, len(dist_list)):
                        temp.append(dist_list[j][t][i])
                    r.append(temp)

    write_c_file(output_file, r)
    

    
    


if __name__ == "__main__": __main__()






        





