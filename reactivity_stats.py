#!/usr/bin/python
#David Tack
#Use a coverage overlap between two (+)DMS libraries to get the overlap you want to keep in the output of this script.
#Won't return stats for any sequence below a length of 10: the deltas here are loopy.

import glob
import numpy
from itertools import islice
import sys

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

def read_in_derived_reactivities(afile):
    '''Reads in a reactivity file'''
    information = {}
    with open(afile) as f:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            information[transcript] = get_reactivity_stats(reactivities)
    return information

def get_reactivity_stats(aline,min_length=10):
    '''Calculates the stats of a line of reactivities'''
    numbers = [float(x) for x in aline.split() if x!= 'NA']
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
    #Generate header
    types,names = ['_max','_average','_std','_gini'],sorted(mega_dict.keys())
    desc = [name+mod for name in names for mod in types]
    header = ['transcript']+desc
    #Combine all keys to a set, get a list to iterate through
    xlist = [z.keys() for z in mega_dict.values()]
    xset = set([j for i in xlist for j in i])
    #Write File
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
    all_dicts = {}
    for fyle in sorted(glob.glob('*.react')):
         temp_dict = read_in_derived_reactivities(fyle)
         all_dicts[fyle.split('.')[0]] = temp_dict
    #
    overlap_keys = get_good_keys(sys.argv[1])
    dump_large_file(all_dicts,overlap_keys,sys.argv[2])


