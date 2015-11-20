#! /usr/bin/env python

def guess_encoding(fh,nsamp=10000):
  '''
    S - Sanger        Phred+33,  raw reads typically (0, 40)
    X - Solexa        Solexa+64, raw reads typically (-5, 40)
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
        with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
  '''
  if isinstance(fh,basestring):
    fh = open(fh,'rU')
  lines = (l.strip('\n') for l in fh)
  minq = 'z'
  maxq = '!'
  for i,l in enumerate(lines):
    if i % 4 == 3:
      minq = min(minq+l)
      maxq = max(maxq+l)
    if i >= nsamp*4:
      break
  if maxq > 'K': # encoding is S or L
    return 'Phred+64'
  else:
    return 'Phred+33'

if __name__ == '__main__':
  import argparse
  import sys
  parser = argparse.ArgumentParser(description='Guess the encoding of fastq file')
  parser.add_argument('infile', nargs='?', type=argparse.FileType('rU'), default=sys.stdin)
  args = parser.parse_args()
  print >>sys.stdout, guess_encoding(args.infile)
