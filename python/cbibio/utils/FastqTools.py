#! /usr/bin/env python
import os

def split_paired_unpaired(infile,pairs,unpaired=None):
  from Fastq import FastqFile
  from collections import OrderedDict
  
  # Setup output files
  outh = open(pairs,'w') if isinstance(pairs,basestring) else pairs
  if unpaired is not None:
    outh2 = open(unpaired,'w') if isinstance(unpaired,basestring) else unpaired
  
  buffer = OrderedDict()
  for s in FastqFile(infile).records():
    if s.basename() not in buffer:
      buffer[s.basename()] = s
    else:
      print >>outh, buffer.pop(s.basename())
      print >>outh, s
      if unpaired is not None:
        for v in buffer.values():
          print >>outh2, v
      buffer = OrderedDict()

def guess_encoding(infile,nsamp=10000):
  '''
    S - Sanger        Phred+33,  raw reads typically (0, 40)
    X - Solexa        Solexa+64, raw reads typically (-5, 40)
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
        with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
  '''
  from Fastq import FastqFile
  qrange = ('z','!')
  for s in FastqFile(infile).records():
    qrange = (min(s.qual+qrange[0]), max(s.qual+qrange[1]))
    nsamp -= 1
    if nsamp <= 0: break
  
  if qrange[1] > 'K': # encoding is S or L
    return 'Phred+64'
  else:
    return 'Phred+33'
