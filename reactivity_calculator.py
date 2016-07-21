#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO
import math
import os
import argparse
from header.__init__ import *
from argparse import RawTextHelpFormatter




        
#Derive raw(before final normalization) reactivity       
def derive_raw_react(trans, distri_p, distri_m):
    result = {}
    for i in range(0, len(trans)):      
        for j in range(len(distri_p[trans[i]])):
            distri_p[trans[i]][j] = math.log((int(distri_p[trans[i]][j])+1),math.e) #take ln to the read counts
        for j in range(len(distri_m[trans[i]])):
            distri_m[trans[i]][j] = math.log((int(distri_m[trans[i]][j])+1),math.e)       
        s_p = sum(distri_p[trans[i]])
        s_m = sum(distri_m[trans[i]])
        length = len(distri_p[trans[i]])
        if s_p!= 0 and s_m!= 0:
            r = []
            for j in range(0, len(distri_p[trans[i]])):
                f_p = (float(distri_p[trans[i]][j]))/float(s_p)*length #normalize by the average logged read counts on each transcript
                f_m = (float(distri_m[trans[i]][j]))/float(s_m)*length
                raw_react = f_p-f_m
                r.append(max(0, raw_react)) #Set negative reactivity as 0
            result[trans[i]] = r
    return result



#Calculate 2-8% normalization scale for each transcript
def cal_norm_level(d_react, transcripts, nt_s):
    result = []
    for t in d_react:
        r = d_react[t]
        a = []
        for k in range(1,len(r)):
            if transcripts[t][k-1] in nt_s:
                a.append(r[k])
        a.sort(reverse = True)
        eight = a[int(len(a)*0.02):int(len(a)*0.1)] #get the rest of top 8% data after removing the top 2%
        meight = sum(eight)/len(eight) #Calculate the average of data (2-8% normalization scale)
        if meight > 0:
            result.append([t, meight])
    return result


#Calculate final reactivity
def cal_final_reactivity(d_react, transcripts, nt_s, threshold, norm_s):
    result = {}
    nls = {}
    for i in range(len(norm_s)):
        nls[norm_s[i][0]] = norm_s[i][1] #Get transcript normalization scale

    threshold = float(threshold)
    for t in d_react:
        r = d_react[t]
        re = []
        if t in nls:
            meight = float(nls[t])
            for k in range(1, len(r)):
                if transcripts[t][k-1] in nt_s:
                    re.append(str(float('%.3f'%min((r[k]/meight), threshold)))) #Raw reactivity on each nucleotide devide by the normalization scale and capped by the specified threshold 
                else:
                    re.append('NA')
            re.append('NA')
            result[t] = re
    return result
                


def write_file(file_name, a):
    h = file(file_name, 'w')
    for i in range(len(a)):
        for j in range(len(a[i])-1):
            h.write(str(a[i][j])+"\t")
        j = j+1
        h.write(str(a[i][j])+"\n")
    h.close()



def __main__():
    parser = argparse.ArgumentParser(usage = "python react_cal_structurefold.py <mode> <RTSC_minus> <RTSC_plus> <Reference_file> <norm_scale> <reactivity_file> [-nt <Specificity>] [-t Thres] [-h]",
                                     description='Derive structural reactivity from RT stop count file\n'+
                                 "Mode 1: Output structural reactivity and 2-8% normalization scale for each transcript\n"+
                                 "Mode 2: Output structural reactivity with input 2-8% normalization scale for each transcript\n", formatter_class=RawTextHelpFormatter)
    parser.add_argument('mo', metavar="<mode>",help='Mode of the program (1 or 2)')
    parser.add_argument('dfm', metavar="<RTSC_minus>",help='(-)DMS/SHAPE RT stop count file')
    parser.add_argument('dfp', metavar="<RTSC_plus>",help='(+)DMS/SHAPE RT stop count file')
    parser.add_argument('seq', metavar="<Reference_file>",help='(Transcriptome) reference library used to map the reads (fasta format)')
    parser.add_argument('-nt', "--nucl", metavar="<Specificity>",default = 'AC', help='Specificity of nucleotides to derive reactivity (AC[default] for DMS, ATCG for SHAPE)')
    parser.add_argument('norm', nargs='?', metavar="<norm_scale>",help='normalization scale for each transcript (output for mode 1; input for mode 2 [Default: [RTSC name]_normalization_scale.txt (Default value only for mode 1)]')
    parser.add_argument('react', nargs='?', metavar="<reactivity_file>",help='Output reactivity file (text file) [Default: [RTSC name].react]')

    parser.add_argument("-t", "--thres", dest="Thres", default = 7, type = float, help="Threshold to cap the reactivity (Default 7)", metavar="Thres")

    args = parser.parse_args()
    #print(args)
    

    
    dist_file1 = args.dfp #plus library
    dist_file2 = args.dfm #minus library
    seq_file = args.seq #Reference library(genome/cDNA)
    nt_spec = args.nucl #only show reactivity for AC or ATCG
    thres = args.Thres #Threshold to cap the reactivities
    output_file = args.react 
    norm_level = args.norm
    mode = args.mo

    dist_p = parse_dist(dist_file1)
    dist_m = parse_dist(dist_file2)

    n_m = dist_file2.split(".rtsc")[0].strip()
    p_m = dist_file1.split(".rtsc")[0].strip()
    prefix = n_m+"_"+p_m

    if output_file == None:
        if mode == '1':
            output_file = prefix+".react"
        else:
            output_file = prefix+"_norm.react"

    if norm_level == None:
        if mode == '1':
            norm_level = prefix+"_normalization_scale.txt"
        else:
            sys.stderr.write("Error: Missing normalization scales from input!\n")
            sys.exit()

        
        

    #syspathrs = os.getcwd()
    h = file(output_file,'w')

    seqs = SeqIO.parse(open(seq_file),'fasta')
    nt = set()
    for i in range(len(nt_spec)):
        nt.add(nt_spec[i])

    flag = 0
    tra = []
    dist_p = dist_p[1]
    dist_m = dist_m[1]
    transcript = {}
    for seq in seqs:
        n = seq.id
        tra.append(n)
        transcript[n] = str(seq.seq)

    derived_reactivity = derive_raw_react(tra, dist_p, dist_m)
    if mode == '1':
        nl = cal_norm_level(derived_reactivity, transcript, nt)
        write_file(norm_level, nl)

    else:
        nl = read_t_file(norm_level)


    fsn = cal_final_reactivity(derived_reactivity, transcript, nt, thres, nl)

    for t in fsn:
        h.write(t+"\n")
        for j in range(len(fsn[t])-1):
            h.write(fsn[t][j]+"\t")
        j = j+1
        h.write(fsn[t][j]+"\n")


    h.close()
        

if __name__ == "__main__" : __main__()


            
    
    
        





















        





