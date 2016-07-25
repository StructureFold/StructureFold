#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from header.read_file import *
from header.write_file import *
from argparse import RawTextHelpFormatter
import argparse
import os
import glob


def combine(lis):
    offset = 0
    r = []
    r.append(lis[0][0][:])
    r.append(lis[0][1][:])
    for i in range(len(lis)):
        oa = int(lis[i][0][0])

        if len(lis[i])>2:
            for j in range(2, len(lis[i])):
                temp = [int(lis[i][j][0])+offset, int(lis[i][j][1])+offset, lis[i][j][2]]
                r.append(temp)
        offset = offset+oa
    return r
                
            

def cal_rmsd(f1, f2):
    f1s = {}
    f2s = {}
    if len(f1)<=2 or len(f2)<=2:
        return 'nan'
    else:
        tps = set()
        s = 0
        c = 0
        for i in range(2, len(f1)):
            tp = (int(f1[i][0]), int(f1[i][1]))
            f1s[tp] = 10**((-1)*float(f1[i][2]))
            tps.add(tp)
        for i in range(2, len(f2)):
            tp = (int(f2[i][0]), int(f2[i][1]))
            f2s[tp] = 10**((-1)*float(f2[i][2]))
            tps.add(tp)

        for tp in tps:
            if tp in f1s:
                f1v = f1s[tp]
            else:
                f1v = 0
            if tp in f2s:
                f2v = f2s[tp]
            else:
                f2v = 0
            s = s+ (f1v-f2v)**2
            c = c+1
        rmsd = (s/c)**0.5
    return rmsd

def parse_direc(p):
    r = p.split('/')[-1].strip()
    if len(r)>0:
        return r
    else:
        return p.split('/')[-2].strip()
            
            
            
            
            
    
        
        
        
    


def __main__():

    parser = argparse.ArgumentParser(description='Calculate RMSD between RNA structures in batch', formatter_class=RawTextHelpFormatter)
    parser.add_argument('id_f', metavar="<RNA_IDs>",help='List of RNAs(IDs) to calculate RMSD')
    parser.add_argument('direc1', metavar="<Directory1>",help='First Directory of output from batch_folder.py (partition function mode) for comparison')
    parser.add_argument('direc2', metavar="<Directory2>",help='Second Directory of output from batch_folder.py (partition function mode) for comparison')
    parser.add_argument('result', nargs='?', metavar="<output_file>",help='Output RMSDs between RNA structures [Default: [Directory1]_vs_[Directory2].rmsd]')


    args = parser.parse_args()

    
    id_list = args.id_f
    direc1 = args.direc1
    direc2 = args.direc2

    if args.result:
        result = args.result
    else:
        result = os.path.basename(direc1)+"_vs_"+os.path.basename(direc2)+".rmsd"
    

    r = []
    ids = read_t_file(id_list)
    for i in range(len(ids)):
        t = ids[i][0]
        flag1 = False
        flag2 = False
        p1 = os.path.join(direc1, "Base-pairing_probability", t+"*")
        p2 = os.path.join(direc2, "Base-pairing_probability", t+"*")
        for f in glob.iglob(p1):
            file1 = f
            flag1 = True
        for f in glob.iglob(p2):
            file2 = f
            flag2 = True
        if flag1 and flag2:
            f1 = read_t_file(file1)
            f2 = read_t_file(file2)
            rmsd = cal_rmsd(f1, f2)
            r.append([t, rmsd])
    write_c_file(result, r)

    

    

    
    


if __name__ == "__main__": __main__()






        





