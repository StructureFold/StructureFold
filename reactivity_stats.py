#!/usr/bin/python
#David Tack
#Returns max, average, standard deviation, and gini of reactivity for each transcript in .react files.
#Operates on the entire directory at once, each set of columns will be named according to .react file names.
#Minumum length to operate on (flag m) is set to 10bp
#Length to ignore from 3' end is set to 0 (flag n)
#Note that these two interact; 


#Imports
import argparse
import glob
import numpy
from itertools import islice

def get_good_keys(afile):
    '''Read in an overlap file to prune for sequences with good coverage'''
    good = {}
    with open(afile,'r') as f:
        for line in f:
            good[line.strip()] = 'NULL'
    return good

def gini(list_of_values):
    '''Returns the GINI value of a list of values'''
    try:
        sorted_list = sorted(list_of_values)
        height, area = 0, 0
        for value in sorted_list:
            height += value
            area += height - value / 2.
        fair_area = height * len(list_of_values) / 2.
        return (fair_area - area) / fair_area
    except ZeroDivisionError:
        return 'NA'

def read_in_derived_reactivities(afile,min_len,trim_numb):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile,'r') as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = get_reactivity_stats(reactivities,min_len,trim_numb)
    return information

def get_reactivity_stats(aline,min_length,trim_number):
    '''Calculates the stats of a line of reactivities'''
    numbers = [float(x) for x in aline.split() if x!= 'NA']
    if trim_number:
        numbers = numbers[:-trim_number]
    ###
    if len(numbers) < min_length:
        numbers_2 = 'NA','NA','NA','NA'
        return numbers_2
    try:
      numbers_2 = max(numbers),numpy.average(numbers),numpy.std(numbers),gini(numbers)
    except ValueError:
      numbers_2 = 'NA','NA','NA','NA'
    return numbers_2


def dump_large_file(mega_dict,keepers,outfile='reactivity_stats.csv'):
    '''Generate File'''
    types,names = ['_max','_average','_std','_gini'],sorted(mega_dict.keys())
    desc = [name+mod for name in names for mod in types]
    header = ['transcript']+desc
    xlist = [z.keys() for z in mega_dict.values()]
    xset = set([j for i in xlist for j in i])
    with open(outfile,'w') as g:
        g.write(','.join(header)+'\n')
        for transcript in xset:
            if transcript in keepers:
                zline = [mega_dict[name][transcript] if transcript in mega_dict[name] else ['NA','NA','NA','NA'] for name in names]
                zlineII = [q for w in zline for q in w]
                zlineIII = [str(p) for p in zlineII]
                g.write(','.join([transcript]+zlineIII)+'\n')

if __name__ == '__main__':
    print ''
    print '\033[1;4;94mStructure Fold:\033[0;0;92m reactivity_stats.py\033[0m'
    print ''
    parser = argparse.ArgumentParser(description='Generates a simple statistical report for all .react files in the directory.')
    parser.add_argument("accepted", type=str, help="list of transcripts to keep from all .react files")
    parser.add_argument('-n',type=int, default=0, help='<int> [default = 20] ignore n last bp of reactivity',dest='trim')
    parser.add_argument('-m',type=int, default=10, help='<int> [default = 10] minumum length of transcript',dest='minlen')
    parser.add_argument('-o',type=str, default='stats_out', help='<str> [default = stats_out] output file ',dest='out')
    args = parser.parse_args()
    all_dicts = {}
    for fyle in sorted(glob.glob('*.react')):
         temp_dict = read_in_derived_reactivities(fyle,args.minlen,args.trim)
         all_dicts[fyle.split('.')[0]] = temp_dict
    overlap_keys = get_good_keys(args.accepted)
    dump_large_file(all_dicts,overlap_keys,args.out+'_'+str(args.trim)+'n_'+str(args.minlen)+'m'+'.csv')


