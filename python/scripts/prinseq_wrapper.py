#! /usr/bin/env python
import sys
from collections import defaultdict
from glob import glob
import re

from cbibio.utils.guess_encoding import guess_encoding


""" Constants """
# Extensions for fastq files
FASTQ_EXTENSIONS = set(['fastq','fq'])
# Extensions for fasta files
FASTA_EXTENSIONS = set(['fasta','fa', 'fas', 'fna'])
# Arguments for PRINSEQ
PRINSEQ_ARGS = ['min_len', 'max_len', 'range_len',
                'min_gc','max_gc','range_gc',
                'min_qual_score', 'max_qual_score',
                'min_qual_mean', 'max_qual_mean',
                'ns_max_p', 'ns_max_n', 
                'seq_num', 'derep', 'derep_min',
                'lc_method', 'lc_threshold',
                'trim_to_len',
                'trim_left', 'trim_right', 
                'trim_left_p','trim_right_p',
                'trim_tail_left','trim_tail_right',
                'trim_ns_left','trim_ns_right',
                'trim_qual_left','trim_qual_right',
                'trim_qual_type', 'trim_qual_rule', 'trim_qual_window', 'trim_qual_step',
                 ]
# Flags for PRINSEQ                 
PRINSEQ_FLAGS = ['noniupac']
# Preset for Illumina data
ILLUMINA_PRESET = { 'min_len':40,
                    'min_qual_mean':20,
                    'derep': '23',
                    'lc_method': 'dust',
                    'lc_threshold':7,
                    'trim_qual_left': 20,
                    'trim_qual_right': 20,
                    'trim_tail_left': 5,
                    'trim_tail_right': 5,
                    'trim_ns_left': 1,
                    'trim_ns_right': 1,
                    }
# Preset for IonTorrent data
                    



""" Functions """
def chunks(l, n):
  ''' Yield successive n-sized chunks from list. '''
  for i in xrange(0, len(l), n):
    yield l[i:i+n]

"""
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
"""

"""
def getcmdlist(allfiles, args):
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
"""

def make_proccmd(fs, args):
    fastx = 'fasta' if args.fasta else 'fastq'
    paired = len(fs) == 2
    
    # Get the prefix of the files
    extstr = '|'.join(FASTA_EXTENSIONS) if args.fasta else '|'.join(FASTQ_EXTENSIONS)    
    m = re.search('^(.*)%s([a-zA-Z0-9]+).(%s)$' % (args.delim, extstr),  fs[0])
    if m:
        prefix = m.group(1)
    else:
        prefix = '.'.join(fs[0].split('.')[:-1])
    
    # Get encoding of the files
    if guess_encoding(fs[0]) == 'Phred+64':
        commandline.append('-phred64')
        
    # Build command line
    commandline = ['prinseq-lite',]    
    commandline.append('-%s' % fastx)
    commandline.append('%s' % fs[0])
    if paired:
        commandline.append('-%s2' % fastx)
        commandline.append('%s' % fs[1])

    if not args.out_prinseq_names:
        commandline.append('-out_good')
        commandline.append('%s_pseq' % prefix)
    else:
        pass # do prinseq default

    if not args.out_bad:
        commandline.append('-out_bad')
        commandline.append('null')
    else:
        if not args.out_prinseq_names:
            commandline.append('-out_bad')
            commandline.append('%s_pseqfail' % prefix)
        else:
           pass # do prinseq default

    commandline.append('$basecmd')
    commandline.append('&')    
    return commandline

def make_basecmd(args):
    ''' Convert wrapper arguments to PRINSEQ command line '''
    commandline = []
    if not args.qual_header: commandline.append('-no_qual_header')
    
    d = vars(args)    
    for argkey in PRINSEQ_ARGS:
        if d[argkey] is not None:
            commandline.append('-%s' % argkey)
            commandline.append('%s' % d[argkey])
        elif args.illumina and argkey in ILLUMINA_PRESET:
            commandline.append('-%s' % argkey)
            commandline.append('%s' % ILLUMINA_PRESET[argkey])
    for argkey in PRINSEQ_FLAGS:
        if d[argkey]:
            commandline.append('-%s' % argkey)        

    return commandline

def get_filesets(flist, args):
    ''' Attempt to pair files 
        Assumes that paired files have the same prefix ahead of some delimiter.
    '''
    if args.single:
        return [(f,) for f in flist]
    pairs = defaultdict(list)
    for f in flist:
        extstr = '|'.join(FASTA_EXTENSIONS) if args.fasta else '|'.join(FASTQ_EXTENSIONS)
        m = re.search('^(.*)%s([a-zA-Z0-9]+).(%s)$' % (args.delim, extstr),  f)
        if m:
            pairs[m.group(1)].append(f)
        else:
            pairs[f].append(f)
    ret = []
    for prefix, fl in pairs.iteritems():
        if len(fl) == 1:
            ret.append((fl[0],))
        elif len(fl) == 2:
            ret.append(tuple(sorted(fl)))
        else:
            sys.exit("ERROR: identifying pairs failed! %s" % ', '.join(fl))
    return ret

def make_header(args):
  ''' Returns header for batch scripts '''
  return [
          '#! /bin/bash',
          '#SBATCH -t %d' % args.walltime,
          '#SBATCH -p %s' % args.partition,
          '#SBATCH -N %d' % args.nodes,
          '',
          'module load %s' % args.prinseq_module,
          'export basecmd="%s"' % ' '.join(make_basecmd(args)),'',
          'echo "[-- prinseq_wrapper --]"',
          'echo "[-- prinseq_wrapper --] Prinseq executable: $(which prinseq-lite)"',
          'echo "[-- prinseq_wrapper --] Prinseq arguments:  $basecmd"','',          
         ]

def make_footer():
  return ['',
          '# Wait for all subprocesses to finish',
          'wait',
          'echo "[-- prinseq_wrapper --] Complete"',
         ]

def main(args):
    if not args.nosubmit:
        from subprocess import Popen,PIPE

    if args.fofn is not None:
        allfiles  = sorted([l.strip() for l in args.fofn])
    else:
        # Generate commands for each file set
        if args.fasta:
            allfiles  = sorted([f for f in glob('*') if f.split('.')[-1] in FASTA_EXTENSIONS])
        else:
            allfiles  = sorted([f for f in glob('*') if f.split('.')[-1] in FASTQ_EXTENSIONS])
        
    filesets = get_filesets(allfiles, args)
    print >>sys.stderr, "INPUT FILES FOUND:"
    for fs in filesets:
        if len(fs)==1:
            print >>sys.stderr, 'unpaired: %s' % fs[0]
        elif len(fs)==2:
            print >>sys.stderr, 'reads1: %s\treads2: %s' % (fs)

    header = make_header(args)
    footer = make_footer()
    
    for i,chunk in enumerate(chunks(filesets, args.chunksize)):
        body = [' '.join(make_proccmd(fs, args)) for fs in chunk]
        script = header + body + footer
        if args.nosubmit:
            with open('job.%02d.sh' % (i+1), 'w') as outh: 
                print >>outh, '\n'.join(script)
        else:
            #submit these
            print >>sys.stderr, '\n### Job %d ###' % (i+1)
            print >>sys.stderr, '\n'.join(script)


"""



    cmdlist = getcmdlist(allfiles, args)

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
"""

if __name__ == '__main__':
  import argparse
  import sys
  parser = argparse.ArgumentParser(description='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  preset_group = parser.add_argument_group("Presets")
  preset_group.add_argument('--illumina', action='store_true',
                           help='''Illumina filtering presets. 
                                   Equivalent to "%s"''' % ' '.join('-%s %s' % (k,v) for k,v in ILLUMINA_PRESET.iteritems()) )
  
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
  slurm_group.add_argument('--prinseq_module', default="prinseq/0.20.4", 
                           help='Name of prinseq module. Calling "module load [prinseq_module]" must load prinseq-lite into environment')
  
  # preset_group = parser.add_argument_group("Presets options")
  
  output_group = parser.add_argument_group("Output options")
  output_group.add_argument('--out_format', type=int, default=3,
                           help='''To change the output format, use one of the following options. 
                                   If not defined, the output format will be the same as the input format. 
                                   1 (FASTA only), 2 (FASTA and QUAL), 3 (FASTQ), 4 (FASTQ and FASTA), or 5 (FASTQ, FASTA and QUAL)''')
  output_group.add_argument('--out_prinseq_names', action='store_true',
                           help='''By default, this wrapper names the good output files by removing the file
                                   extension and adding "_pseq" to the file name. The prinseq default
                                   is to add random characters to prevent overwriting; however, the
                                   filenames can be cumbersome for downstream analysis. Use this flag to
                                   use the default prinseq behavior.''')
  """                                 
  output_group.add_argument('--out_good', 
                           help='''By default, this wrapper names the good output files by removing the file
                                   extension and adding "_prinseq" to the file name. The prinseq default
                                   is Note that this is 
                                   different the prinseq default, which 
                                   By default, the output files are created in the same directory
                                   as the input file containing the sequence data with an
                                   additional "_prinseq_good_XXXX" in their name (where XXXX is
                                   replaced by random characters to prevent overwriting previous
                                   files). To change the output filename and location, specify the
                                   filename using this option. The file extension will be added
                                   automatically (either .fasta, .qual, or .fastq). For paired-end
                                   data, filenames contain additionally "_1", "_1_singletons",
                                   "_2", and "_2_singletons" before the file extension. Use
                                   "-out_good null" to prevent the program from generating the
                                   output file(s) for data passing all filters. Use "-out_good
                                   stdout" to write data passing all filters to STDOUT (only for
                                   FASTA or FASTQ output files). Example: use "file_passed" to generate the output file
                                   file_passed.fasta in the current directory)''')
  """
  output_group.add_argument('--out_bad', action='store_true',
                            help='''By default, this wrapper does not output "bad" data that
                                    does not pass filters. Use this flag to output the failed
                                    data; the file names will have "_pseqfail" appended to the file name,
                                    unlessed used in conjunction with "--out_prinseq_names"''')
                                    
  output_group.add_argument('--qual_header', action='store_true',
                            help='''By default, this wrapper outputs an empty header line for the
                                    quality data to reduce file size. Use this flag to enable the
                                    quality header.''')

  filter_group = parser.add_argument_group("Filter options")
  filter_group.add_argument('--min_len', type=int,
                            help=" Filter sequence shorter than min_len.")
  filter_group.add_argument('--max_len', type=int,
                            help=" Filter sequence longer than max_len.")
  filter_group.add_argument('--range_len', 
                            help='''Filter sequence by length range. Multiple range values
                                    should be separated by comma without spaces. Example:
                                    --range_len 50-100,250-300''')
  filter_group.add_argument('--min_gc',
                            help='''Filter sequence with GC content below min_gc.''')
  filter_group.add_argument('--max_gc',
                            help='''Filter sequence with GC content above max_gc.''')
  filter_group.add_argument('--range_gc', 
                            help='''Filter sequence by GC content range. Multiple range values
                                    should be separated by comma without spaces. Example:
                                    --range_gc 50-100,250-300''')

  filter_group.add_argument('--min_qual_score', 
                            help='''Filter sequence with at least one quality score below min_qual_score.''')                                  
  filter_group.add_argument('--max_qual_score', 
                            help='''Filter sequence with at least one quality score above max_qual_score.''')
  filter_group.add_argument('--min_qual_mean', 
                            help='''Filter sequence with quality score mean below min_qual_mean.''')
  filter_group.add_argument('--max_qual_mean', 
                            help='''Filter sequence with quality score mean above max_qual_mean.''')                            

  filter_group.add_argument('--ns_max_p', 
                            help='''Filter sequence with more than ns_max_p percentage of Ns.''')                                       
  filter_group.add_argument('--ns_max_n', 
                            help='''Filter sequence with more than ns_max_n Ns.''')

  filter_group.add_argument('--noniupac', 
                            help='''Filter sequence with characters other than A, C, G, T or N.''')
  filter_group.add_argument('--seq_num', 
                            help='''Only keep the first seq_num number of sequences (that pass all other filters).''')

  filter_group.add_argument('--derep',
                            help='''Type of duplicates to filter. Allowed values are 1, 2,
                                    3, 4 and 5. Use integers for multiple selections (e.g.
                                    124 to use type 1, 2 and 4). The order does not
                                    matter. Option 2 and 3 will set 1 and option 5 will
                                    set 4 as these are subsets of the other option. 1
                                    (exact duplicate), 2 (5' duplicate), 3 (3' duplicate),
                                    4 (reverse complement exact duplicate), 5 (reverse
                                    complement 5'/3' duplicate)''')
  filter_group.add_argument('--derep_min',
                            help='''This option specifies the number of allowed
                                    duplicates. If you want to remove sequence duplicates
                                    that occur more than x times, then you would specify
                                    x+1 as the -derep_min values. For examples, to remove
                                    sequences that occur more than 5 times, you would
                                    specify -derep_min 6. This option can only be used in
                                    combination with -derep 1 and/or 4 (forward and/or
                                    reverse exact duplicates).''')
  filter_group.add_argument('--lc_method',
                            help='''Method to filter low complexity sequences. The current
                                    options are "dust" and "entropy". Use "-lc_method
                                    dust" to calculate the complexity using the dust
                                    method''')
  filter_group.add_argument('--lc_threshold',
                            help='''The threshold value (between 0 and 100) used to filter
                                    sequences by sequence complexity. The dust method uses
                                    this as maximum allowed score and the entropy method
                                    as minimum allowed value.''')

  trim_group = parser.add_argument_group("Trim options")
  trim_group.add_argument('--trim_to_len', type=int,
                          help='''Trim all sequence from the 3'-end to result in
                                  sequence with this length.''')
  trim_group.add_argument('--trim_left', type=int,
                          help='''Trim sequence at the 5'-end by trim_left positions.''')
  trim_group.add_argument('--trim_right', type=int,
                          help='''Trim sequence at the 3'-end by trim_right positions.''')
  trim_group.add_argument('--trim_left_p', type=int,
                          help='''Trim sequence at the 5'-end by trim_left_p percentage
                                  of read length. The trim length is rounded towards the
                                  lower integer (e.g. 143.6 is rounded to 143
                                  positions). Use an integer between 1 and 100 for the
                                  percentage value.''')
  trim_group.add_argument('--trim_right_p', type=int,
                          help='''Trim sequence at the 3'-end by trim_right_p percentage
                                  of read length. The trim length is rounded towards the
                                  lower integer (e.g. 143.6 is rounded to 143
                                  positions). Use an integer between 1 and 100 for the
                                  percentage value.''')
  trim_group.add_argument('--trim_tail_left', type=int,
                          help='''Trim poly-A/T tail with a minimum length of
                                  trim_tail_left at the 5'-end.''')
  trim_group.add_argument('--trim_tail_right', type=int,
                          help='''Trim poly-A/T tail with a minimum length of
                                  trim_tail_right at the 3'-end.''')
  trim_group.add_argument('--trim_ns_left', type=int,
                          help='''Trim poly-N tail with a minimum length of trim_ns_left
                                  at the 5'-end.''')
  trim_group.add_argument('--trim_ns_right', type=int,
                          help='''Trim poly-N tail with a minimum length of
                                  trim_ns_right at the 3'-end.''')
  trim_group.add_argument('--trim_qual_left', type=int,
                          help='''Trim sequence by quality score from the 5'-end with
                                  this threshold score.''')
  trim_group.add_argument('--trim_qual_right', type=int,
                          help='''Trim sequence by quality score from the 3'-end with
                                  this threshold score.''')
  trim_group.add_argument('--trim_qual_type', choices=['min', 'mean', 'max', 'sum'],
                          help='''Type of quality score calculation to use. Allowed
                                  options are min, mean, max and sum.''')
  trim_group.add_argument('--trim_qual_rule', choices=['lt','gt','et'],
                          help='''Rule to use to compare quality score to calculated
                                  value. Allowed options are lt (less than), gt (greater
                                  than) and et (equal to).''')
  trim_group.add_argument('--trim_qual_window', type=int,
                          help='''The sliding window size used to calculate quality
                                  score by type. To stop at the first base that fails
                                  the rule defined, use a window size of 1.''')
  trim_group.add_argument('--trim_qual_step', type=int,
                          help='''Step size used to move the sliding window. To move the
                                  window over all quality scores without missing any,
                                  the step size should be less or equal to the window
                                  size.''')

  args = parser.parse_args()
  main(args)


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

ILLUMINA_PARAMS = {
    
}
