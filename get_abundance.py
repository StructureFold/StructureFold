#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from header.parse_dis_react import *
from argparse import RawTextHelpFormatter
import argparse




def __main__():
    parser = argparse.ArgumentParser(description='Calculate abundance of each transcript (normalized by length)\n', formatter_class=RawTextHelpFormatter)
    parser.add_argument('rtsc', metavar="<RTSC_file>",help='input RT stop count file')
    parser.add_argument('result', nargs='?', metavar="<abundance_file>",help='Output abundance file (text file) [Default: [RTSC name]_abundance.txt]')
    args = parser.parse_args()
    
    dist_file1 = args.rtsc
    output_file = args.result


    n_m = dist_file1.split(".rtsc")[0].strip()
    prefix = n_m+"_"

    if output_file == None:
        output_file = prefix+"abundance.txt"

        
    dist1 = parse_dist(dist_file1)
    dist1 = dist1[1]
    dis1_sum = {}
    for t in dist1:
        dis1_sum[t] = 0
        for i in range(len(dist1[t])):
            dis1_sum[t] = dis1_sum[t] + int(dist1[t][i])
    with open(output_file, 'w') as h:
        for t in dis1_sum:
            h.write(t+"\t"+str(float(dis1_sum[t])/len(dist1[t])*1000)+"\n")
    
    
    
    


if __name__ == "__main__": __main__()






        





