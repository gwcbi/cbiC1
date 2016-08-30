#! /usr/bin/env python
import sys
from collections import defaultdict
from glob import glob
import re

from cbibio.utils.guess_encoding import guess_encoding
from cbibio.utils.findfiles import sizeof_fmt

""" Constants """

# Arguments for PRINSEQ
PRINSEQ_ARGS =        [ 'min_len', 'max_len', 'range_len',
                        'min_gc','max_gc','range_gc',
                        'min_qual_score', 'max_qual_score',
                        'min_qual_mean', 'max_qual_mean',
                        'ns_max_p', 'ns_max_n', 
                        'seq_num',
                        'derep', 'derep_min',
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
ILLUMINA_PRESET =     { 'min_len': 40,
                        'min_qual_mean': 20,
                        'derep': '23',
                        'lc_method': 'dust',
                        'lc_threshold': 7,
                        'trim_qual_left': 20,
                        'trim_qual_right': 20,
                        'trim_tail_left': 5,
                        'trim_tail_right': 5,
                        'trim_ns_left': 1,
                        'trim_ns_right': 1,
                      }
# Preset for IonTorrent data
IONTORRENT_PRESET =   { 'min_len': 40,
                        'min_qual_mean': 17,
                        'derep': '14',
                        'lc_method': 'dust',
                        'lc_threshold': 7,
                        'trim_qual_left': 17,
                        'trim_qual_right': 17,
                        'trim_tail_left': 5,
                        'trim_tail_right': 5,
                        'trim_ns_left': 1,
                        'trim_ns_right': 1,
                      }

""" Functions """

def chunks(l, n):
  ''' Yield successive n-sized chunks from list. '''
  for i in xrange(0, len(l), n):
    yield l[i:i+n]

def get_file_prefix(fn, args):
    extstr = '|'.join(FASTA_EXTENSIONS) if args.fasta else '|'.join(FASTQ_EXTENSIONS)    
    m = re.search('^(.*)%s([a-zA-Z0-9]+)\.(%s)$' % ('_', extstr),  fn)
    if m:
        return m.group(1)
    m = re.search('^(.*)\.(%s)$' % extstr, fn)
    if m:
        return m.group(1)
    return '.'.join(fn.split('.')[:-1])

def make_filecmd(rfset, args):
    ''' Creates the command line for ReadFileSet '''
    encoding = '-phred64' if rfset.encoding() == 'Phred+64' else ''
    if rfset.unpaired():
        if rfset.UP.zipped:
            file1str = '-%s <(gunzip -c %s)' % (rfset.UP.filetype.lower(), rfset.UP.fullpath)
        else:
            file1str = '-%s %s' % (rfset.UP.filetype.lower(), rfset.UP.fullpath)
        file2str = ''
    elif rfset.paired() or rfset.both():
        # Read 1
        if rfset.R1.zipped:
            file1str = '-%s <(gunzip -c %s)' % (rfset.R1.filetype.lower(), rfset.R1.fullpath)
        else:
            file1str = '-%s %s' % (rfset.R1.filetype.lower(), rfset.R1.fullpath)
        # Read 2
        if rfset.R2.zipped:
            file2str = '-%s2 <(gunzip -c %s)' % (rfset.R2.filetype.lower(), rfset.R2.fullpath)
        else:
            file2str = '-%s2 %s' % (rfset.R2.filetype.lower(), rfset.R2.fullpath)
    if rfset.both():
        # Concatenate unpaired reads to read1 file
        if rfset.R1.zipped:
            file1str = '-%s <(cat %s %s | gunzip -c )' % (rfset.R1.filetype.lower(), rfset.R1.fullpath, rfset.UP.fullpath)
        else:
            file1str = '-%s <(cat %s %s)' % (rfset.R1.filetype.lower(), rfset.R1.fullpath, rfset.UP.fullpath)
    
    prefix = rfset.UP.basename if rfset.unpaired else rfset.R1.basename
    # Output arguments
    if not args.out_prinseq_names:
        goodstr = "-out_good %s.pseq" % prefix             # Output good file with wrapper naming scheme
    else:
        goodstr = ''                                       # Output good file with prinseq naming scheme
    
    if not args.out_bad:
        badstr = '-out_bad null'                           # Do not output bad file
    else:
        if not args.out_prinseq_names:
            badstr =  '-out_bad %s.pseqfail' % prefix      # Output bad file with wrapper naming scheme
        else:
            badstr =  ''                                   # Output bad file with prinseq naming scheme
    commandline = ['prinseq-lite', encoding, file1str, file2str, goodstr, badstr, "$basecmd", ]
    return commandline

def make_basecmd(args):
    ''' Convert wrapper arguments to PRINSEQ command line '''
    commandline = []
    if not args.qual_header: commandline.append('-no_qual_header')
    
    d = vars(args)    
    for argkey in PRINSEQ_ARGS:
        if d[argkey] is not None:
            commandline.append('-%s %s' % (argkey, d[argkey]))
        elif args.illumina and argkey in ILLUMINA_PRESET:
            commandline.append('-%s %s' % (argkey, ILLUMINA_PRESET[argkey]))        
        elif args.iontorrent and argkey in IONTORRENT_PRESET:
            commandline.append('-%s %s' % (argkey, IONTORRENT_PRESET[argkey]))        
    for argkey in PRINSEQ_FLAGS:
        if d[argkey]:
            commandline.append('-%s' % argkey)
    return commandline

def make_header(args):
  ''' Returns header for batch scripts '''
  return [     '#! /bin/bash',
               '#SBATCH -t %d' % args.walltime,
               '#SBATCH -p %s' % args.partition,
               '#SBATCH -N 1',
               '',
               'SN="prinseq_wrapper"',
               'echo "[---$SN---] ($(date)) Starting $SN"',
               't1=$(date +"%s")',
               '',
               'module load %s' % args.prinseq_module,
               'export basecmd="%s"' % ' '.join(make_basecmd(args)),'',
               'echo "[---$SN---] ($(date))  Prinseq executable: $(which prinseq-lite)"',
               'echo "[---$SN---] ($(date))  Prinseq arguments:  $basecmd"',
               '',          
         ]

def make_footer():
  ''' Returns footer for batch scripts '''
  return [     '',
               '# Wait for all subprocesses to finish',
               'wait',
               '',
               '#---Complete job',
               't2=$(date +"%s")',
               'diff=$(($t2-$t1))',
               'echo "[---$SN---] Total time: ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds."',
               'echo "[---$SN---] ($(date)) $SN COMPLETE."',
         ]

from cbibio.utils.findfiles import BioFile, ReadFileSet, rec_find_files, filter_files
from cbibio.utils.findfiles import find_mate_files_byreadname, find_mate_files_byfilename
from cbibio.utils.chunkfiles import distribute_readsets
import os

def main(args):
    # Get the list of files
    if args.fofn is not None:
        fiter = (l.strip() for l in args.fofn)
    elif args.infiles == '.':
        fiter = rec_find_files('.')
    else:
        fiter = args.infiles

    biofiles = [BioFile(f) for f in fiter if os.path.isfile(f)]
    biofiles = filter_files(biofiles, 'FASTA' if args.fasta else 'FASTQ', args.include_files, args.exclude_files)
    if len(biofiles) == 0: sys.exit("No matching files were found.")
    print >>sys.stderr, '%d %s files were found' % (len(biofiles), 'FASTA' if args.fasta else 'FASTQ')

    if args.single:
        print >>sys.stderr, 'Assuming all files are unpaired...'
        rfsets = [ReadFileSet(up=bf) for bf in biofiles]
    else:
        print >>sys.stderr, 'Finding mate files...'
        if args.by_readname:
            rfsets = find_mate_files_byreadname(biofiles)
        else:
            rfsets = find_mate_files_byfilename(biofiles, not args.no_check)

    print >>sys.stderr, '\nReadFile sets:'
    print >>sys.stderr, '\n'.join(str(rfset) for rfset in rfsets)    
    print >>sys.stderr, '%d file sets were found:' % (len(rfsets))
    print >>sys.stderr, '\t%d paired sets' % sum(rfset.paired() for rfset in rfsets)
    print >>sys.stderr, '\t%d unpaired sets' % sum(rfset.unpaired() for rfset in rfsets)
    print >>sys.stderr, '\t%d paired+unpaired sets' % sum(rfset.both() for rfset in rfsets)  

    if not args.quiet:
        response = None
        while response is None:
            response = raw_input("Do you wish to continue? [Y/N] ")
            if response in set(['N','n','no','NO','No']):
                sys.exit("Aborted")
            elif response not in set(['Y','y','yes','Yes','YES']):
                response = None
        assert response in set(['Y','y','yes','Yes','YES'])

    nodes = distribute_readsets(rfsets, args.nslots, args.nsequential, args.njobs)
    
    # Construct job script
    header = make_header(args)
    footer = make_footer()
    for i,node in enumerate(nodes):
        print >>sys.stderr, '#################################### Node%02d ####################################' % i
        s = '## nproc: %d size: %s' % (len(node), sizeof_fmt(sum(_.size() for _ in node)))
        print >>sys.stderr, '%s##' % s.ljust(78)
        for rfset in node:
            for l in str.split('\n'):
                s = '## %s' % l
                print >>sys.stderr, '%s##' % s.ljust(78)
        print >>sys.stderr, '################################################################################' 
        
        
        body = ['module load parallel',
                'parallel -j %d <<EOF' % args.nslots,
               ]
        for rfset in node:
            body.append( ' '.join(make_filecmd(rfset, args)) )
        body.append('EOF')
        script = header + body + footer
        if args.nosubmit:
            with open('job.%02d.sh' % (i+1), 'w') as outh: 
                print >>outh, '\n'.join(script)
            print >>sys.stderr, '[--- job: prinseq%02d ---] Slurm script written to "job.%02d.sh".' % ((i+1), (i+1))
        else:
            # Submit
            from subprocess import Popen,PIPE
            p = Popen(['sbatch','-J','prinseq%02d' % (i+1)],  stdin=PIPE, stdout=PIPE)
            out,err = p.communicate(input='\n'.join(script))
            print >>sys.stderr, '[--- job: prinseq%02d ---] %s. Running %d prinseq processes.' % ((i+1), out.strip('\n'), len(node))
        if args.detail:
            print >>sys.stderr, '\n'.join(script)
            print >>sys.stderr, ''

if __name__ == '__main__':
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--quiet', action='store_true')
    
    preset_group = parser.add_argument_group("Presets")
    preset_group.add_argument('--illumina', action='store_true',
                             help='''Illumina filtering presets. 
                                     Equivalent to "%s"''' % ' '.join('-%s %s' % (k,v) for k,v in ILLUMINA_PRESET.iteritems()) )
    preset_group.add_argument('--iontorrent', action='store_true',
                             help='''IonTorrent filtering presets. 
                                     Equivalent to "%s"''' % ' '.join('-%s %s' % (k,v) for k,v in IONTORRENT_PRESET.iteritems()) )
    
    input_group = parser.add_argument_group("Input options")
    input_group.add_argument('--fofn', type=argparse.FileType('r'),
                             help="File of file names")
    input_group.add_argument('--single', action='store_true',
                             help="All files are single-end. Do not attempt to find paired files.")
    input_group.add_argument('--fasta', action='store_true',
                             help="Input files are fasta format.")
    input_group.add_argument('--include_files',
                             help='''Only include files that match the given search string. Can
                                     be a regular expression.''')
    input_group.add_argument('--exclude_files',
                             help='''Exclude files that match the given search string. Can be a
                                     regular expression.''')
    input_group.add_argument('--by_readname', action='store_true',
                              help='''Use readname instead of filename to find pairs. The
                                      advantage here is that the files can have any random name,
                                      as long as the read names match they will be paired. There
                                      are some pitfalls that might cause this to fail: 1) If you
                                      have multiple copies of the same file with different
                                      filenames; and 2) If some file sets include both paired and 
                                      unpaired files.''')
    input_group.add_argument('--no_check', action='store_true',
                              help='''Do not check the readnames when finding pairs by 
                                      filename. Ignored if --by_readname provided.''')
    
    slurm_group = parser.add_argument_group("Slurm options")
    slurm_group.add_argument('--nosubmit', action='store_true',
                             help="Do not submit jobs directly. Output is written to batch scripts")
    slurm_group.add_argument('--walltime', type=int, default=480, 
                             help="Slurm walltime request (in minutes)")
    slurm_group.add_argument('--partition', default="short", 
                             help="Slurm partition request for resource allocation")
    # slurm_group.add_argument('--nodes', type=int, default=1, 
    #                          help='''Slurm request for nodes per job.''')
    slurm_group.add_argument('--njobs', type=int,
                             help='''Total number of SLURM jobs to start. Each job will use
                                     one node.''')
    slurm_group.add_argument('--nslots', type=int, default=16, 
                             help='''Number of simultaneous processes that can run on a
                                     node. For single-threaded programs, this is equal to the
                                     number of CPUs on each node. Otherwise, nslots should be
                                     the total number of CPUs divided by the number of CPUs
                                     used by each process (rounded up). For example, if each
                                     process uses 300% CPU and each node has 16 CPUs, nslots
                                     should be set to 5''')
    slurm_group.add_argument('--nsequential', type=int, default=1, 
                             help='''Number of process to run sequentially on the same node.
                                     If queue wait times are long, it may be advantageous to
                                     run sequentially to avoid waiting for additional
                                     allocations. For example, if you have 32 read sets to
                                     process, and --processes_per_node=16, the wrapper will
                                     request two allocations. If --sequential_processes=2, 
                                     then only one allocation will be requested.''')
    slurm_group.add_argument('--prinseq_module', default="prinseq/0.20.4", 
                             help='Name of prinseq module. Calling "module load [prinseq_module]" must load prinseq-lite into environment')
    slurm_group.add_argument('--detail', action='store_true',
                             help='Show detailed information about jobs being submitted.')
    
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument('--out_format', type=int, default=3,
                             help='''To change the output format, use one of the following options. 
                                     If not defined, the output format will be the same as the input format. 
                                     1 (FASTA only), 2 (FASTA and QUAL), 3 (FASTQ), 4 (FASTQ and FASTA), or 5 (FASTQ, FASTA and QUAL)''')
    output_group.add_argument('--out_prinseq_names', action='store_true',
                             help='''By default, this wrapper names the good output files by removing the file
                                     extension and adding ".pseq" to the file name. The prinseq default
                                     is to add random characters to prevent overwriting; however, the
                                     filenames can be cumbersome for downstream analysis. Use this flag to
                                     use the default prinseq behavior.''')
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
    
    parser.add_argument('infiles', nargs='*', default=".", 
                        help='''Input files. File name expansion of wildcards is handled by
                                the OS. Recursively searches current directory if no 
                                argument is present.''')
    args = parser.parse_args()
    main(args)
