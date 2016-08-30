#! /usr/bin/env python

# This is the set of characters that are ONLY found in phred33
# !"#$%&\'()*+,-./0123456789:
phred33 = set(chr(_) for _ in range(33,59)) 



def guess_encoding(fh,nsamp=10000):
    '''
    S - Sanger        Phred+33,  raw reads typically (0, 40).  ASCII 33-73
    X - Solexa        Solexa+64, raw reads typically (-5, 40). ASCII 59-104
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40).  ASCII 64-104
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40).  ASCII 66-104
        with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold) 
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)   ASCII 33-74
    '''
    if isinstance(fh,basestring):
        fh = open(fh,'rU')
    lines = (l for l in fh)
    minq = 'z'
    maxq = '!'
    for i,l in enumerate(lines):
        if i % 4 == 3:
            # Solexa+64 can be as low as -5, so if there are any ASCII characters below 
            # 59 (64-5), it is definitely phred-33.
            # chr(64-5) == ';'
            if any(c < ';' for c in l.strip('\n')):
                return 'Phred+33'
            # Phred+33 maxes out at 74 (33+41) for Illumina data, but PacBio QV can be
            # above 60. We'll assume that we will not see any QVs above 64. If there are 
            # ASCII characters above 97 (33+64) we will assume this is Phred+64.
            # chr(33+64) == 'a'
            if any(c >= 'a' for c in l.strip('\n')):
                return 'Phred+64'
            # Otherwise just set the overall min and max values
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
