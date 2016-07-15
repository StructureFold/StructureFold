#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from header.parse_dis_react import *
from header.write_file import *
from argparse import RawTextHelpFormatter
import argparse




def __main__():
    parser = argparse.ArgumentParser(description='Calculate correlation between biological replicate libraries\n', formatter_class=RawTextHelpFormatter)
    parser.add_argument('rep1', metavar="<RTSC_file1>",help='input RT stop count file of replicate 1')
    parser.add_argument('rep2', metavar="<RTSC_file2>",help='input RT stop count file of replicate 2')
    parser.add_argument('result', nargs='?', metavar="<correlation_file>",help='Output abundance file (text file) [Default: [RTSC name]_correlation.txt]')
    args = parser.parse_args()
    
    dist_file1 = args.rep1
    dist_file2 = args.rep2
    output_file = args.result


    n_m1 = dist_file1.split(".rtsc")[0].strip()
    n_m2 = dist_file2.split(".rtsc")[0].strip()
    prefix = n_m1+"_"+n_m2

    if output_file == None:
        output_file = prefix+"_correlation.txt"

        
    dist1 = parse_dist(dist_file1)
    dist1 = dist1[1]

    dist2 = parse_dist(dist_file2)
    dist2 = dist2[1]

    r = []
    for t in dist1:
        for i in range(len(dist1[t])):
            temp = [t, dist1[t][i], dist2[t][i]]
            r.append(temp)

    write_t_file(output_file, r)
    
    
    
    


if __name__ == "__main__": __main__()






        





