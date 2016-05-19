#! /usr/bin/env python

import sys
from collections import defaultdict
from glob import glob
import re
import gzip

from cbibio.utils.guess_encoding import guess_encoding

""" Constants """
# Extensions for fastq files
FASTQ_EXTENSIONS = set(['fastq','fq'])
# Extensions for fasta files
FASTA_EXTENSIONS = set(['fasta','fa', 'fas', 'fna'])

def get_file_prefix(fn, args):
    # Strip off gzip
    prefix = re.sub('\.gz$', '', fn, flags=re.I)
    extensions = FASTA_EXTENSIONS if args.fasta else FASTQ_EXTENSIONS
    for ext in extensions:
        prefix = re.sub('\.%s$' % ext, '', prefix, flags=re.I)
    
    if not args.single:
        if args.delim in prefix:
            prefix = args.delim.join(prefix.split(args.delim)[:-1])
    return prefix

ILLUMINA_PRESET =     { 'trim_qual': 0.05,
                        'drop_ambig': True,
                        'min_len': 40,
                        'fastx_artifacts': True,
                        'fastx_qual': 20,
                        'fastx_qual_pct': 50,
                      }

SCRIPT_HEADER = '''#! /bin/bash
#SBATCH -t %(walltime)s
#SBATCH -p %(partition)s
#SBATCH -N %(nodes)s

SN="seqtk-wrapper"
echo "[---$SN---] ($(date)) Starting $SN"
t1=$(date +"%%s")

module load cbiC1
module load %(seqtk_module)s
module load %(fastx_module)s
echo "[---$SN---] ($(date))  seqtk executable: $(which seqtk)"

'''

SCRIPT_FOOTER = '''
#--- Wait for all subprocesses to finish
wait

#---Complete job
t2=$(date +"%%s")
diff=$(($t2-$t1))
echo "[---$SN---] Total time: ($(date)) $(($diff / 60)) minutes and $(($diff %% 60)) seconds."
echo "[---$SN---] ($(date)) $SN COMPLETE."

'''

""" Functions """
def chunks(l, n):
  ''' Yield successive n-sized chunks from list. '''
  for i in xrange(0, len(l), n):
    yield l[i:i+n]

def get_filesets(flist, args):
    ''' Attempt to pair files 
        Assumes that paired files have the same prefix ahead of some delimiter.
    '''
    if args.single:
        return [(f,) for f in flist]
    pairs = defaultdict(list)
    for f in flist:
        prefix = get_file_prefix(f, args)
        pairs[prefix].append(f)
    
    ret = []
    for prefix, fl in pairs.iteritems():
        if len(fl) == 1:
            ret.append((fl[0],))
        elif len(fl) == 2:
            ret.append(tuple(sorted(fl)))
        else:
            sys.exit("ERROR: identifying pairs failed! %s" % ', '.join(fl))
    return sorted(ret, key=lambda x:x[0])

def make_filecmd(fs, args):
    pipeline = []

    gzipped = fs[0].endswith('.gz')
    paired = len(fs) == 2
    prefix = get_file_prefix(fs[0], args)

    # Start input stream using cat (for SE) or mergepe (for PE)
    if paired:
        pipeline.append('seqtk mergepe %s %s' % (fs[0],fs[1]))
    else:
        if gzipped:
            pipeline.append('gunzip -c %s' % fs[0])
        else:
            pipeline.append('cat %s' % fs[0])

    # Adjust the encoding of the files
    if gzipped:
        encoding = guess_encoding(gzip.open(fs[0]))
    else:
        encoding = guess_encoding(fs[0])
    # print >>sys.stderr, '%s %s' % (fs[0],encoding)
    if encoding == 'Phred+64':
        pipeline.append('seqtk seq -Q64 -V -')
    
    # End trimming of reads
    if args.trim_left or args.trim_right:
        nleft  = args.trim_left if args.trim_left else 0
        nright = args.trim_right if args.trim_right else 0        
        pipeline.append('seqtk trimfq -b %d -e %d -' % (nleft, nright))
    
    # Quality trimming of reads
    if args.trim_qual:
        pipeline.append('seqtk trimfq -q %f -' % (args.trim_qual))

    # Filtering
    argC = '' if args.keep_comment else '-C'
    argL = '-L %d' % args.min_len if args.min_len else ''
    argN = '-N' if args.drop_ambig else ''
    pipeline.append(' '.join(['seqtk', 'seq', argC, argL, argN, '-']))
    
    # Fastq quality filter
    if args.fastx_qual:
        pipeline.append('fastq_quality_filter -Q33 -q %d -p %d' % (args.fastx_qual, args.fastx_qual_pct))

    # Fastq artifact filter
    if args.fastx_artifacts:
        pipeline.append('fastx_artifacts_filter -Q33')
    
    # Output
    if paired:
        out1 = '%s.seqtk%s1.fq' % (prefix, args.delim)
        out2 = '%s.seqtk%s2.fq' % (prefix, args.delim)
        pipeline.append('seqtk dropse')
        pipeline.append('deinterleave_fastq %s %s' % (out1,out2))
    else:
        out1 = '%s.seqtk.fq' % (prefix)
        pipeline[-1] = '%s > %s' % (pipeline[-1], out1)

    return ' | '.join(pipeline)

def main(args):
    if args.fofn is None:
        allfiles = []    
        extensions = FASTA_EXTENSIONS if args.fasta else FASTQ_EXTENSIONS
        for ext in extensions:
            allfiles.extend(glob('*.%s' % ext))
            allfiles.extend(glob('*.%s.gz' % ext))
    else:
        allfiles  = sorted([l.strip() for l in args.fofn])

    # Identify pairs of files
    filesets = get_filesets(allfiles, args)
    print >>sys.stderr, "INPUT FILES FOUND:"
    for fs in filesets:
        if len(fs)==1:
            print >>sys.stderr, 'unpaired: %s' % fs[0]
        elif len(fs)==2:
            print >>sys.stderr, 'reads1: %s\treads2: %s' % (fs)
    print >>sys.stderr, ''

    for i,chunk in enumerate(chunks(filesets, args.chunksize)):
        filecmds = []
        for fs in chunk:
            cmd = '%s &' % make_filecmd(fs, args)
            filecmds.extend(['echo "[---$SN---] ($(date)) COMMAND: %s"' % cmd, cmd, ''])
        
        script = [SCRIPT_HEADER % vars(args)] + filecmds + [SCRIPT_FOOTER % vars(args)]
        if args.nosubmit:
            with open('job.%02d.sh' % (i+1), 'w') as outh: 
                print >>outh, '\n'.join(script)
            print >>sys.stderr, '[--- job: seqtk%02d ---] Slurm script written to "job.%02d.sh".' % ((i+1), (i+1))                
        else:
            from subprocess import Popen,PIPE
            p = Popen(['sbatch','-J','seqtk%02d' % (i+1)],  stdin=PIPE, stdout=PIPE)
            out,err = p.communicate(input='\n'.join(script))
            # out = 'Dummy slurm response 12345.\n'            
            print >>sys.stderr, '[--- job: seqtk%02d ---] %s. Running %d seqtk processes.' % ((i+1), out.strip('\n'), len(chunk))            
        
        if args.detail:
            print >>sys.stderr, '\n'.join(script)
            print >>sys.stderr, ''  

if __name__ == '__main__':
  import argparse
  import sys
  parser = argparse.ArgumentParser(description='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  preset_group = parser.add_argument_group("Presets")
  preset_group.add_argument('--illumina', action='store_true',
                            help='''Illumina filtering presets. 
                                    Equivalent to "%s"''' % ' '.join('--%s %s' % (k,v) for k,v in ILLUMINA_PRESET.iteritems()) )
  # preset_group.add_argument('--iontorrent', action='store_true',
  #                          help='''IonTorrent filtering presets. 
  #                                  Equivalent to "%s"''' % ' '.join('-%s %s' % (k,v) for k,v in IONTORRENT_PRESET.iteritems()) )
  
  input_group = parser.add_argument_group("Input options")
  input_group.add_argument('--fofn', type=argparse.FileType('r'),
                           help="File of file names")
  input_group.add_argument('--single', action='store_true',
                           help="All files are single-end. Do not attempt to find paired files.")
  input_group.add_argument('--delim', default='_',
                           help="Delimiter used when looking for paired files.")
  input_group.add_argument('--fasta', action='store_true',
                           help="Input files are fasta format.")
  
  slurm_group = parser.add_argument_group("Slurm options")
  slurm_group.add_argument('--nosubmit', action='store_true',
                           help="Do not submit jobs directly. Output is written to batch scripts")
  slurm_group.add_argument('--chunksize', type=int, default=16, 
                           help="Number of prinseq instances per job. This should be less than or equal to the number of CPUs per node")
  slurm_group.add_argument('--walltime', type=int, default=480, 
                           help="Slurm walltime request (in minutes)")
  slurm_group.add_argument('--partition', default="short", 
                           help="Slurm partition request for resource allocation")
  slurm_group.add_argument('--nodes', type=int, default=1, 
                           help="Slurm node request")
  slurm_group.add_argument('--seqtk_module', default="seqtk/1.0-r68-git", 
                           help='Name of seqtk module. Calling "module load [seqtk_module]" must load seqtk into environment')
  slurm_group.add_argument('--fastx_module', default="fastx/0.0.13", 
                           help='Name of fastx module. Calling "module load [fastx_module]" must load fastx binaries into environment')
  slurm_group.add_argument('--detail', action='store_true',
                           help='Show detailed information about jobs being submitted.')
  
  output_group = parser.add_argument_group("Output options")
  output_group.add_argument('--keep_comment', action='store_true',
                           help='''Keep comments at the header lines''')
  
  filter_group = parser.add_argument_group("Filter options")
  filter_group.add_argument('--min_len', type=int,
                            help=" Filter sequence shorter than min_len.")
  filter_group.add_argument('--drop_ambig', action='store_true',
                            help=" Drop sequences containing ambiguous bases.")
  filter_group.add_argument('--min_qual_score', 
                            help='''Filter sequence with at least one quality score below min_qual_score.''')                                  
  
  trim_group = parser.add_argument_group("Trim options")
  trim_group.add_argument('--trim_left', type=int,
                          help='''Trim sequence at the 5'-end by trim_left positions.''')
  trim_group.add_argument('--trim_right', type=int,
                          help='''Trim sequence at the 3'-end by trim_right positions.''')
  trim_group.add_argument('--trim_qual', type=float,
                          help='''Trim using phred algorithm with error rate threshold.''')

  other_group = parser.add_argument_group("Other options")
  other_group.add_argument('--fastx_artifacts', action='store_true',
                          help='''Use FASTX toolkit to filter artifacts.''')
  other_group.add_argument('--fastx_qual', type=int,
                          help='''Use FASTX toolkit quality filter. fastx_qual is the cutoff for determining low-quality bases.''')                          
  other_group.add_argument('--fastx_qual_pct', type=int,
                          help='''Use FASTX toolkit quality filter. fastx_qual_pct is the minimum percent of bases 
                                  that must have low quality for read to be dropped. Only activated if fastx_qual is set.''')

  args = parser.parse_args()

  # Deal with presets
  if args.illumina:
      for k,v in ILLUMINA_PRESET.iteritems():
          if getattr(args,k) is None: # Argument was not set
              setattr(args, k, v)
          elif getattr(args,k) is False: # Argument default is False
              setattr(args, k, v)

  main(args)
