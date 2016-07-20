#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from header.parse_dis_react import *
from argparse import RawTextHelpFormatter
import argparse




def __main__():
    parser = argparse.ArgumentParser(description='Combine RT stops from biological replicate libraries\n', formatter_class=RawTextHelpFormatter)
    parser.add_argument('rep1', metavar="<RTSC_file1>",help='input RT stop count file of replicate 1')
    parser.add_argument('rep2', metavar="<RTSC_file2>",help='input RT stop count file of replicate 2')
    parser.add_argument('result', nargs='?', metavar="<RTSC_file>",help='Output RTSC file (text file) [Default: [RTSC name]_combine.txt]')
    args = parser.parse_args()
    
    dist_file1 = args.rep1
    dist_file2 = args.rep2
    output_file = args.result


    n_m1 = dist_file1.split(".rtsc")[0].strip()
    n_m2 = dist_file2.split(".rtsc")[0].strip()
    prefix = n_m1+"_"+n_m2

    if output_file == None:
        output_file = prefix+"_combine.txt"

        
    dist1 = parse_dist(dist_file1)
    dist1 = dist1[1]

    dist2 = parse_dist(dist_file2)
    dist2 = dist2[1]

    r = {}
    for t in dist1:
        if t not in r:
            r[t] = []
        for i in range(len(dist1[t])):
            r[t].append(int(dist1[t][i])+int(dist2[t][i]))

    with open(output_file, 'w') as h:
        for t in r:
            h.write(t+"\n")
            for i in range(len(r[t])-1):
                h.write(str(r[t][i])+"\t")
            i = i+1
            h.write(str(r[t][i])+"\n")
            h.write("\n")
            
        
    
    
    
    


if __name__ == "__main__": __main__()






        





