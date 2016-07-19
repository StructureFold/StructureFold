#!/usr/bin/python
#David Tack
#Part of a set of tools to simply use of Structure Fold
#Run in a folder with all of your coverage files. Generates flat text files of the overlap 
#between each pair of coverage files, beyond a certain threshold value, default 1

#Imports
import glob
import os
import argparse

#Functions
def get_coverage(afyle):
    '''Reads in a coverage file'''
    xdict = {}
    with open(afyle,'r') as f:
        for line in f:
            name, value = line.strip().split('\t')
            xdict[name]= float(value)
    return xdict

def filter_nested_dictionaries(adict,threshold):
    '''Removes dictionary entries with values less than a given threshold'''
    for sub_dict in adict.values():
        for k, v in sub_dict.items():
            if v <= threshold :
                del sub_dict[k]

def inspect_setwise(adict):
    '''Generate all possible overlaps'''
    dump_dict={}
    for k in sorted(adict.keys()):
        sub_list = list(adict.keys())
        print '####',k
        for item in sorted(sub_list):
            aself = set(adict[k].keys())
            other = set(adict[item].keys())
            print item, len(aself.intersection(other))
            dump_dict[(k,item)] = aself.intersection(other)
        print ''
    return dump_dict

def set_dict_to_files(adict,threshold):
    '''Create a directory, write all overlap files in that directory'''
    os.mkdir('Overlaps_'+str(threshold))
    os.chdir('Overlaps_'+str(threshold))
    for k, v in adict.items():
        with open('x'.join([k[0],k[1]])+'.txt','w') as g:
            for item in v:
                g.write(item+'\n')

if __name__ == '__main__':
    print ''
    print '\033[1;4;94mStructure Fold:\033[0;0;92m generate_coverages.py\033[0m'
    print ''
    parser = argparse.ArgumentParser(description='Creates transcript lists with coverage between data sets.')
    parser.add_argument('-n',type=int, default=1, help='<int> [default = 1] coverage threshold to use',dest='threshold')
    args = parser.parse_args()
    all_dicts= {}
    for fyle in sorted(glob.glob('*_coverage.txt')):
        temp_dict = get_coverage(fyle)
        all_dicts[fyle.split('_')[0]] = temp_dict
    filter_nested_dictionaries(all_dicts,args.threshold)
    transcript_sets = inspect_setwise(all_dicts)
    set_dict_to_files(transcript_sets,args.threshold)
