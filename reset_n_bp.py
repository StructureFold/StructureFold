#!/usr/bin/python
#David Tack
#resets the last n bp of react files in the directory. 
#useful for when the last n bases have low confidence or are unresolvable by the nature of the technique used

#Imports
import glob
from itertools import islice
import argparse

#Function
def reset_last_n_bp(afile,trim_numb):
    '''Writes a new file with the last n bp reset to NA'''
    newfyle = afile.split('.')[0]+'_'+str(trim_numb)+'trim.react'
    with open(afile,'r') as f,open(newfyle,'w') as g:
        while True:
            next_n_lines = list(islice(f, 2))
            if not next_n_lines:
                break
            transcript,reactivities = [n.strip() for n in next_n_lines]
            reacts = [x for x in reactivities.split()]
            if len(reacts) <= trim_numb:
                out_reacts = ['NA']*len(reacts)
            else:
                base_reacts = reacts[:-trim_numb]
                mod_reacts = ['NA']*trim_numb
                out_reacts = base_reacts+mod_reacts
            g.write(transcript+'\n')
            g.write('\t'.join(out_reacts)+'\n')

if __name__ == '__main__':
    print ''
    print '\033[1;4;94mStructure Fold:\033[0;0;92m reset_n_bp.py\033[0m'
    print ''
    parser = argparse.ArgumentParser(description='sets the last n bp of all .react files in the directory to NA')
    parser.add_argument('-n',type=int, default=20, help='<int> [default = 20] last n bp to be reset to NA',dest='trim')
    args = parser.parse_args()
    for fyle in sorted(glob.glob('*.react')):
      reset_last_n_bp(fyle,args.trim)
