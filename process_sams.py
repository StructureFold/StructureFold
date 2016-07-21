#!/usr/bin/env python
import sys
import glob
import subprocess
import time

#Classes
#TBH, you should make a formal class and actually deploy it to both deal with first mis-matches, and the orientation of the read
class Mapped_Read(object):
  def __init__(self,qname='',flag='',rname='',pos='',mapq='',cigar='',rnext='',pnext='',tlen='',seq='',qual='',attributes=''):
    '''Mapped Illumina Read'''
    self.qname = qname
    self.flag = flag
    self.rname = rname
    self.pos = pos
    self.mapq = mapq
    self.cigar = cigar
    self.rnext = rnext
    self.pnext = pnext
    self.tlen = tlen
    self.seq = seq
    self.qual = qual
    self.attributes = attributes

  def __str__(self):
    return ' '.join ([self.name,str(self.paired),self.bitflag,self.cigar,str(self.start),str(self.mapped),str(self.gaps)])

  @classmethod
  def from_sam_line(klass,sam_line):
    '''SAM line -> Mapped_Read'''
    bits= sam_line.strip().split()
    if bits[0].startswith('@'):
      return klass('')
    QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL = bits[0:11]
    if len(bits)==11:
      return klass(QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL)
    if len(bits) > 11:
      ATTR = bits[11:]
      return klass(QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL,ATTR)

def derp_function(infile):
  '''Docstring'''
  with open(infile) as f:
    for line in f:
      a = Mapped_Read.from_sam_line(line)
      #print a.attributes
      print a.cigar, len(a.seq), a.attributes
      #You dummy, it works like a CIGAR string, where the first flag is just numbers
      #If header line,write
      #If line with mismatch at 1
      #if line where bitflag aligns to the opposite stand

#Functions
def remove_first_mismatch(infile):
  '''docstring'''
  outfile= infile.split('.')[0]+'_NoSNP0.sam'
  with open(infile) as f,open(outfile,'w') as g:
    for line in f:
      flag,flag2,tl=0,0,line.strip().split('\t')
      for i in range(len(tl)):
        if tl[i].strip().find("NM:i:")!=-1:
          n,flag2 = tl[i].strip().split(':')[-1].strip(),1
          if int(n) > 3:
            flag = 1
        if tl[i].strip().find("MD:Z:")!=-1:
          md,flag2 = tl[i].strip().split(':')[-1].strip(),1
          if (int(md[0])) == 0:
            flag = 1
      if flag2 == 0:
        continue
      if flag == 0:
        g.write(line)
  return outfile

def keep_good_reads(infile):
  '''docstring'''
  outfile= infile.split('.')[0]+'_Fwd.sam'
  #samtools view -h -F 0xeff -S AA_cdna.sam > AA_2.sam
  subprocess.call(' '.join(['samtools','view','-h','-F','0xeff','-S',infile,'>',outfile]),shell=True)
  return outfile

#Workflow
if __name__ == '__main__':
  print ''
  print '\033[1;4;94mRemoving Improper Reads and Mismatches at postion 0 for all .SAM in directory\033[0m'
  print ''
  start_time,count = time.asctime(),0
  for fyle in sorted(glob.glob('*.sam')):
    print '\033[92mProcessing with SAMtools:\033[0m {}'.format(fyle)
    out1 = keep_good_reads(fyle)
    print '\033[92mManually removing SNPs in position 0 from:\033[0m {}'.format(out1)
    out2 = remove_first_mismatch(out1)
    print '\033[92mWrote final output:\033[0m {}'.format(out2)
    print '\033[91mRemoving Intermediate file:\033[0m {}'.format(out1)
    subprocess.call(' '.join(['rm',out1]),shell=True)
    print ''
    count+=1
  print '\033[1;4;94mStart Time:\033[0m ',start_time
  print '\033[1;4;94mEnd Time:\033[0m ',time.asctime()
  print '\033[1;4;92mTotal Files Processed:\033[0m ',count
  print ''
    
    
    
