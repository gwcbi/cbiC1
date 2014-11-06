#! /usr/bin/env python
'''
Usage:   seqtk trimfq [options] <in.fq>

Options: -q FLOAT    error rate threshold (disabled by -b/-e) [0.05]
         -l INT      maximally trim down to INT bp (disabled by -b/-e) [30]
         -b INT      trim INT bp from left (non-zero to disable -q/-l) [0]
         -e INT      trim INT bp from right (non-zero to disable -q/-l) [0]
'''


'''
Usage:   seqtk seq [options] <in.fq>|<in.fa>

Options: -q INT    mask bases with quality lower than INT [0]
         -X INT    mask bases with quality higher than INT [255]
         -n CHAR   masked bases converted to CHAR; 0 for lowercase [0]
         -l INT    number of residues per line; 0 for 2^32-1 [0]
         -Q INT    quality shift: ASCII-INT gives base quality [33]
         -s INT    random seed (effective with -f) [11]
         -f FLOAT  sample FLOAT fraction of sequences [1]
         -M FILE   mask regions in BED or name list FILE [null]
         -L INT    drop sequences with length shorter than INT [0]
         -c        mask complement region (effective with -M)
         -r        reverse complement
         -A        force FASTA output (discard quality)
         -C        drop comments at the header lines
         -N        drop sequences containing ambiguous bases
         -1        output the 2n-1 reads only
         -2        output the 2n reads only
         -V        shift quality by '(-Q) - 33'
'''

'''
usage: fastq_quality_filter [-h] [-v] [-q N] [-p N] [-z] [-i INFILE] [-o OUTFILE]
Part of FASTX Toolkit 0.0.13 by A. Gordon (gordon@cshl.edu)

   [-h]         = This helpful help screen.
   [-q N]       = Minimum quality score to keep.
   [-p N]       = Minimum percent of bases that must have [-q] quality.
   [-z]         = Compress output with GZIP.
   [-i INFILE]  = FASTA/Q input file. default is STDIN.
   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.
   [-v]         = Verbose - report number of sequences.
                  If [-o] is specified,  report will be printed to STDOUT.
                  If [-o] is not specified (and output goes to STDOUT),
                  report will be printed to STDERR.
'''

'''
usage: fastx_artifacts_filter [-h] [-v] [-z] [-i INFILE] [-o OUTFILE]
Part of FASTX Toolkit 0.0.13 by A. Gordon (gordon@cshl.edu)

   [-h]         = This helpful help screen.
   [-i INFILE]  = FASTA/Q input file. default is STDIN.
   [-o OUTFILE] = FASTA/Q output file. default is STDOUT.
   [-z]         = Compress output with GZIP.
   [-v]         = Verbose - report number of processed reads.
                  If [-o] is specified,  report will be printed to STDOUT.
                  If [-o] is not specified (and output goes to STDOUT),
                  report will be printed to STDERR.
'''
from cbibio.utils.FastqTools import split_paired_unpaired,guess_encoding
from subprocess import Popen,PIPE

def main(parser):
  args = parser.parse_args()
  encoding = guess_encoding(args.fastq)
  paired = args.fastq2 is not None
  
  if paired and args.unpaired and args.out is None:
    print >>sys.stderr, "Error: Output prefix [--out] must be specified for unpaired read output [--unpaired]."
    parser.print_help()
    sys.exit()
  
  if args.verbose:
    print >>sys.stderr, "Encoding: %s" % encoding
    print >>sys.stderr, "Paired:   %s" % paired
    print >>sys.stderr, "Fastq1:   %s" % args.fastq
    if paired:
      print >>sys.stderr, "Fastq2:   %s" % args.fastq2
    if args.out is None:
      print >>sys.stderr, "Output:   STDOUT"
    else:
      print >>sys.stderr,   "Output:   %s.fastq" % args.out
      if args.unpaired:
        print >>sys.stderr, "Unpaired: %s.unpaired.fastq" % args.out
    print >>sys.stderr, '\n'   

  if paired:
    ''' Usage: seqtk mergepe <in1.fq> <in2.fq> '''  
    merge_args = ['seqtk','mergepe', args.fastq, args.fastq2]
    if args.verbose: print >>sys.stderr, ' '.join(merge_args)    
    p0A = Popen(merge_args, stdout=PIPE)
    p0A_out = p0A.stdout
  else:
    p0A = None
    p0A_out = open(args.fastq,'rU')
  
  if encoding == 'Phred+64':
    ''' Usage:   seqtk seq [options] <in.fq>|<in.fa> '''
    enc_args = ['seqtk','seq','-Q64','-V','-']
    if args.verbose: print >>sys.stderr, ' '.join(enc_args)
    p0B = Popen(enc_args, stdin=p0A_out, stdout=PIPE)
    p0A_out.close()
    p0B_out = p0B.stdout
  else:
    p0B = None
    p0B_out = p0A_out    

  ''' Usage:   seqtk trimfq [options] <in.fq> '''
  trim_args = ['seqtk','trimfq','-q',str(args.trim_qual),'-']
  if args.verbose: print >>sys.stderr, ' '.join(trim_args)
  p2 = Popen(trim_args,stdin=p0B_out,stdout=PIPE)
  p0B_out.close()

  ''' Usage:   seqtk seq [options] <in.fq>|<in.fa> '''
  lf_args = ['seqtk','seq','-L', str(args.min_length)]
  if not args.keep_comment: lf_args.append('-C')
  if not args.keep_ambig: lf_args.append('-N')
  lf_args.append('-')
  if args.verbose: print >>sys.stderr, ' '.join(lf_args)
  p3 = Popen(lf_args,stdin=p2.stdout,stdout=PIPE)
  p2.stdout.close()

  ''' usage: fastq_quality_filter [-h] [-v] [-q N] [-p N] [-z] [-i INFILE] [-o OUTFILE] '''
  qf_args = ['fastq_quality_filter','-Q33','-q',str(args.low_qual),'-p',str(args.percent)]
  if args.verbose: print >>sys.stderr, ' '.join(qf_args)
  p4 = Popen(qf_args,stdin=p3.stdout,stdout=PIPE)
  p3.stdout.close()

  ''' usage: fastx_artifacts_filter [-h] [-v] [-z] [-i INFILE] [-o OUTFILE] '''
  art_args = ['fastx_artifacts_filter','-Q33']
  if args.verbose: print >>sys.stderr, ' '.join(art_args)
  p5 = Popen(art_args,stdin=p4.stdout,stdout=PIPE)
  p4.stdout.close()

  if not paired:
    # Nothing left to do except print out the reads
    if args.out is None:
      print >>sys.stdout, p5.communicate()[0]
    else:    
      with open('%s.fastq' % args.out, 'w') as outh:
        print >>outh, p5.communicate()[0]
  else:
    # Unpaired reads are not kept
    if not args.unpaired:
      ''' Usage: seqtk dropSE <in.fq> '''
      drop_args = ['seqtk','dropse']
      if args.verbose: print >>sys.stderr, ' '.join(drop_args)    
      p6 = Popen(drop_args,stdin=p5.stdout,stdout=PIPE)
      p5.stdout.close()
      if args.out is None:
        print >>sys.stdout, p6.communicate()[0]
      else:
        with open('%s.fastq' % args.out, 'w') as outh:
          print >>outh, p6.communicate()[0]
    else:
      split_paired_unpaired(p5.stdout, '%s.fastq' % args.out, unpaired='%s.unpaired.fastq' % args.out)
      p5.stdout.close()
      
if __name__ == '__main__':
  import argparse
  import sys
  parser = argparse.ArgumentParser(description='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  ''' Arguments for seqtk trimfq '''
  parser.add_argument('-v','--verbose',action="store_true",help="Amount of feedback printed to screen")

  ''' Arguments for seqtk trimfq '''
  parser.add_argument('-q','--trim_qual',type=float,default=0.05,help="error rate threshold for Phred trimming algorithm")

  ''' Arguments for seqtk seq '''
  parser.add_argument('-C','--keep_comment',action="store_true", help="do not drop comments at the header lines")
  parser.add_argument('-N','--keep_ambig',action="store_true", help="do not drop sequences containing ambiguous bases")
  parser.add_argument('-L','--min_length',type=int,default=40, help="drop sequences with length shorter than INT")
  
  ''' Arguments for fastq_quality_filter '''
  parser.add_argument('-l','--low_qual', type=int, default=20, help="Cutoff for determining low-quality bases.")
  parser.add_argument('-p','--percent', type=int, default=50, help="Minimum percent of bases that must have [-l] quality.")
  
  ''' No arguments for fastx_artifacts_filter '''

  ''' Output file for unpaired reads '''  
  parser.add_argument('-u','--unpaired', action="store_true", help="Output unpaired reads. File output prefix [--out] must be specified")
  parser.add_argument('-o','--out',help="Prefix for output files")

  parser.add_argument('fastq',help="fastq file name")
  parser.add_argument('fastq2',nargs='?',help="paired fastq file name")

  main(parser)
