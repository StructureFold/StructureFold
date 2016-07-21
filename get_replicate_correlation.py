#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from header.parse_dis_react import *
from header.write_file import *
from argparse import RawTextHelpFormatter
import argparse


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
    args = parser.parse_args()
    
    file_list = args.replicate
    output_file = args.result

    if len(file_list)<2:
        sys.exit("Need at least 2 replicates for correlation!")



    if output_file == None:
        prefix = get_name(file_list)
        output_file = prefix+"_correlation.csv"

    dist_list = []
    for i in range(len(file_list)):
        dist = parse_dist(file_list[i])
        dist = dist[1]
        dist_list.append(dist)



    r = []
    t_set = find_common(dist_list)
    for t in t_set:
        for i in range(len(dist_list[0][t])):
            temp = [t, dist_list[0][t][i]]
            for j in range(1, len(dist_list)):
                temp.append(dist_list[j][t][i])
            r.append(temp)

    write_c_file(output_file, r)
    

    
    


if __name__ == "__main__": __main__()






        





