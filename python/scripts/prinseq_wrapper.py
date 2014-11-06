#! /usr/bin/env python
from glob import glob
from cbibio.utils.guess_encoding import guess_encoding

""" Constants """
# Extensions for fastq files
FASTQ_EXTENSIONS = set(['fastq','fq'])

""" Functions """
def chunks(l, n):
  ''' Yield successive n-sized chunks from list. '''
  for i in xrange(0, len(l), n):
    yield l[i:i+n]

def make_base_cmd(args):
  ''' Make command line with standard parameters '''
  params = ['out_bad','min_len','min_qual_mean','trim_qual_left','trim_qual_right',
            'ns_max_n','derep','trim_tail_left','trim_tail_right','lc_method',
            'lc_threshold']
  d = vars(args)
  pcmd = ' '.join('-%s %s' % (p,d[p]) for p in params)
  # Only have one flag for now
  fcmd = '-no_qual_header' if not args.qual_header else ''
  return '%s %s' % (pcmd,fcmd)

def make_header(args):
  ''' Returns header for batch scripts '''
  return [
          '#! /bin/bash',
          '#SBATCH -t %d:00:00' % args.wallhrs,
          '#SBATCH -p %s' % args.partition,
          '#SBATCH -N %d' % args.nodes,'',
          'export prinseq=%s' % args.prinseq_path,
          'export basecmd="%s"' % make_base_cmd(args),'',
          'echo "[-- prinseq_wrapper --]"',
          'echo "[-- prinseq_wrapper --] Prinseq executable: $prinseq"',
          'echo "[-- prinseq_wrapper --] Prinseq arguments:  $basecmd"','',          
         ]

def make_footer():
  return ['','# Wait for all subprocesses to finish','wait']

def getcmdlist(allfiles,check_pairs,delim):
  ''' Formats command line arguments for list of files
      If check_pairs=True, the filenames (minus extension) are split on the delimiter and
        files with identical prefixes are considered to be pairs.
      The appropriate command line for each set of files is generated, using -fastq2 if
        necessary.
      Encoding offset is also checked and added to command line
  '''
  from collections import defaultdict
  pairs = defaultdict(list)
  for f in allfiles:
    noext = '.'.join(f.split('.')[:-1])
    if check_pairs:
      prefix = delim.join(noext.split(delim)[:-1])
    else:
      prefix = noext
    pairs[prefix].append(f)
  
  cmdlist = []
  for k,files in pairs.iteritems():
    assert len(files)<3, "ERROR: too many files found! %s -> %s" % (k,files)
    cmd = ['$prinseq','-fastq',files[0]]
    if len(files)==2:
      cmd.append('-fastq2')
      cmd.append(files[1])
    if guess_encoding(files[0]) == 'Phred+64':
      cmd.append('-phred64')
    cmd.extend(['-out_good','%s_prinseq ' % k,'$basecmd','&'])
    # cmdlist.append('echo "[-- prinseq_wrapper --] %s"' % ' '.join(cmd))
    cmdlist.append(' '.join(cmd))    
  return cmdlist


def main(parser):
  args = parser.parse_args()
  if not args.nosubmit:
    from subprocess import Popen,PIPE

  # Generate commands for each file set
  allfiles  = sorted([f for f in glob('*') if f.split('.')[-1] in FASTQ_EXTENSIONS])
  check_pairs = not args.single
  cmdlist = getcmdlist(allfiles,check_pairs,args.delim)

  # Generate header and footer lines
  header  = make_header(args)
  footer  = make_footer()
  
  # Produce jobs
  counter = 1
  for chunk in chunks(cmdlist,args.chunksize):
    script = header + chunk + footer
    if args.nosubmit:
      with open('job%02d.sh' % counter,'w') as outh:
        print >>outh, '\n'.join(script)
    else:
      p = Popen(['sbatch','-J','prinseq%02d' % counter],stdin=PIPE,stdout=PIPE)
      out,err = p.communicate(input='\n'.join(script))
      print >>sys.stderr, "%s. Running %d prinseq processes" % (out,len(chunk))
    counter += 1

if __name__ == '__main__':
  import argparse
  import sys
  parser = argparse.ArgumentParser(description='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--fofn',help="File of file names")
  parser.add_argument('--single',action='store_true',help="All files are single-end. Do not attempt to find paired files.")
  parser.add_argument('--nosubmit',action='store_true',help="Do not submit jobs directly. Output is written to batch scripts")
  parser.add_argument('--delim',default='_',help="Delimiter used when looking for paired files.")
  parser.add_argument('--chunksize',type=int,default=16,help="Number of processes per job. This should be equal to the number of CPUs per node")
  
  parser.add_argument('--wallhrs',type=int,default=8,help="Slurm walltime request (hours)")
  parser.add_argument('--partition',default="short",help="Slurm partition request for resource allocation")
  parser.add_argument('--nodes',type=int,default=1,help="Slurm node request")
  parser.add_argument('--prinseq_path',default="/groups/cbi/prinseq-lite-0.20.3/prinseq-lite",help="Path to prinseq executable")
  
  parser.add_argument('--out_bad', default='null', help="File to write data that fails filters.")
  parser.add_argument('--min_len', default='40', help=" Filter sequence shorter than min_len.")
  parser.add_argument('--min_qual_mean', default='30', help="Filter sequence with quality score mean below min_qual_mean.")
  parser.add_argument('--trim_qual_left', default='20', help="Trim sequence by quality score from the 5'-end with this threshold score.")
  parser.add_argument('--trim_qual_right', default='20', help="Trim sequence by quality score from the 3'-end with this threshold score.")
  parser.add_argument('--ns_max_n', default='0', help="Filter sequence with more than ns_max_n Ns.")
  parser.add_argument('--derep', default='23', help="Type of duplicates to filter.")
  parser.add_argument('--trim_tail_left', default='5', help="Trim poly-A/T tail with a minimum length of trim_tail_left at the 5'-end.")
  parser.add_argument('--trim_tail_right', default='5', help="Trim poly-A/T tail with a minimum length of trim_tail_left at the 3'-end.")
  parser.add_argument('--lc_method', default='dust', help="Method to filter low complexity sequences.")
  parser.add_argument('--lc_threshold', default='7', help="The threshold value (between 0 and 100) used to filter sequences by sequence complexity.")
  # Flags
  parser.add_argument('--qual_header',action='store_true')

  main(parser)


"""
PARAMS = {
           'out_bad':'null',           # File to write data that fails filters.
           'min_len':'40',             # Filter sequence shorter than min_len.
           'min_qual_mean':'25',       # Filter sequence with quality score mean below min_qual_mean. 
           'trim_qual_left':'25',      # Trim sequence by quality score from the 5'-end with this threshold score.
           'trim_qual_right':'25',     # Trim sequence by quality score from the 3'-end with this threshold score.
           'ns_max_n':'0',             # Filter sequence with more than ns_max_n Ns.
           'derep':'1',                # Type of duplicates to filter.
           'trim_tail_left':'5',       # Trim poly-A/T tail with a minimum length of trim_tail_left at the 5'-end.
           'trim_tail_right':'5',      # Trim poly-A/T tail with a minimum length of trim_tail_left at the 3'-end.
           'lc_method':'dust',         # Method to filter low complexity sequences.
           'lc_threshold':'7',         # The threshold value (between 0 and 100) used to filter sequences by sequence complexity.
         }
# Default flags for prinseq
FLAGS = [
           'no_qual_header',           # Generate an empty header line for the quality data in FASTQ files
        ] 
"""
"""
def getbasecmd():
  params = ' '.join('-%s %s' % t for t in PARAMS.iteritems())
  flags  = ' '.join('-%s' % t for t in FLAGS)  
  return '$prinseq %s %s ' % (params,flags)
"""
