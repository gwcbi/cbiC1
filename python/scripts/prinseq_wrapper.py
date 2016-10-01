#! /usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Matthew L. Bendall'

import sys
import os

from cbibio.utils import findfiles
from cbibio.utils import chunkfiles

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

SCRIPT_HEADER = '''#! /bin/bash
#SBATCH -t %(walltime)s
#SBATCH -p %(partition)s
#SBATCH -N 1

SN="prinseq_wrapper"
echo "[---$SN---] ($(date)) Starting $SN"
t1=$(date +"%%s")

module load %(prinseq_module)s
echo "[---$SN---] ($(date))  Prinseq executable: $(which prinseq-lite)"

'''

SCRIPT_FOOTER = '''
#--- Wait for all subprocesses to finish
wait

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] Total time: ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds."
echo "[---$SN---] ($(date)) $SN COMPLETE."

'''

""" Subcommands in pipeline for one ReadSet """
PRINSEQ_CMD = '''prinseq-lite %(encoding)s %(goodstr)s %(badstr)s $basecmd '''
UNZIP1 = 'tmp1=$(mktemp) && gunzip -c %(filepath1)s > $tmp1 && '
UNZIP2 = 'tmp2=$(mktemp) && gunzip -c %(filepath2)s > $tmp2 && '
F1FLAG = ' -%(filetype)s %(filepath1)s'
F2FLAG = ' -%(filetype)s2 %(filepath2)s'
F1FLAG_TMP = ' -%(filetype)s $tmp1'
F2FLAG_TMP = ' -%(filetype)s2 $tmp2'
RMF1 = ' && rm $tmp1'
RMF2 = ' && rm $tmp2'
CATU = 'tmp1=$(mktemp) && cat %(filepath1)s %(filepathU)s > $tmp1 && '
UNZIP_CATU = 'tmp1=$(mktemp) && cat %(filepath1)s %(filepathU)s | gunzip -c > $tmp1 && '

""" String together subcommands depending on paired, unpaired, or mixed and zipped or unzipped """
UPAIR_CMD = PRINSEQ_CMD + F1FLAG
UPAIR_ZIP_CMD = UNZIP1 + PRINSEQ_CMD + F1FLAG_TMP + RMF1
PAIR_CMD = PRINSEQ_CMD + F1FLAG + F2FLAG
PAIR_ZIP_CMD = UNZIP1 + UNZIP2 + PRINSEQ_CMD + F1FLAG_TMP + F2FLAG_TMP + RMF1 + RMF2
MIX_CMD = CATU + PRINSEQ_CMD + F1FLAG_TMP + F2FLAG + RMF1
MIX_ZIP_CMD = UNZIP_CATU + UNZIP2 + PRINSEQ_CMD + F1FLAG_TMP + F2FLAG_TMP + RMF1 + RMF2

def fileset_command(rfs, args):
    ''' Creates the command line for ReadFileSet '''
    d = {}
    # Encoding
    d['encoding'] = '-phred64' if rfs.encoding() == 'Phred+64' else ''

    # Output prefix
    indir, bn = rfs.output_prefix()
    if args.outdir:
        prefix = os.path.join(args.outdir, bn)
    else:
        prefix = os.path.join(indir, bn)

    # Output names
    d['goodstr'] = '' if args.out_prinseq_names else "-out_good %s.pseq" % prefix
    if not args.out_bad:
        d['badstr'] = '-out_bad null'                           # Do not output bad file
    else:
        d['badstr'] = '' if args.out_prinseq_names else '-out_bad %s.pseqfail' % prefix
    
    d['filetype'] = rfs.filetype.lower()
    
    if rfs.unpaired:
        d['filepath1'] = rfs.UP.fullpath
        if not rfs.zipped:
            return UPAIR_CMD % d
        else:
            return UPAIR_ZIP_CMD % d        
    elif rfs.paired:
        d['filepath1'] = rfs.R1.fullpath
        d['filepath2'] = rfs.R2.fullpath
        if not rfs.zipped:
            return PAIR_CMD % d
        else:
            return PAIR_ZIP_CMD % d
    elif rfs.mixed:
        d['filepath1'] = rfs.R1.fullpath
        d['filepath2'] = rfs.R2.fullpath
        d['filepathU'] = rfs.U.fullpath
        if not rfs.zipped:
            return MIX_CMD % d
        else:
            return MIX_ZIP_CMD % d

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

def main(args):
    # Get the ReadSets
    if args.fofn:
        infiles = [l.strip('\n') for l in args.fofn if not l.startswith('#')]
    else:
        infiles = args.infiles

    pair_method = 'readname' if args.by_readname else 'filename'
    rfsets = findfiles.find_readsets(infiles, args.file_type, args.include_files,
                                     args.exclude_files, pair_method, not args.no_check, 
                                     args.single, not args.quiet)
    if not args.quiet:
        findfiles.promptabort()
    
    # Assign processes to nodes
    nodes = chunkfiles.distribute_readsets(rfsets, args.nslots, args.nsequential, args.njobs)    
    
    # Build the submit scripts
    header = SCRIPT_HEADER % vars(args)
    footer = SCRIPT_FOOTER
    for i,node in enumerate(nodes):
        if not args.quiet:
            print >>sys.stderr, '#################################### Node%02d ####################################' % (i+1)
            s = '## nproc: %d size: %s' % (len(node), findfiles.fmtsize_h(sum(_.size() for _ in node)))
            print >>sys.stderr, '%s##' % s.ljust(78)
            for rfset in node:
                for l in str.split('\n'):
                    s = '## %s' % l
                    print >>sys.stderr, '%s##' % s.ljust(78)
            print >>sys.stderr, '################################################################################' 

        body = ['export basecmd="%s"' % ' '.join(make_basecmd(args)),
                'module load parallel',
                'parallel -j %d <<EOF' % args.nslots
        ]
        for rfs in node:
            body.append(fileset_command(rfs, args))
        
        body.append('EOF')
        job_sh = [header] + body + [footer]
        
        if args.nosubmit:
            with open('prinseq_job.%02d.sh' % (i+1), 'w') as outh: 
                print >>outh, '\n'.join(job_sh)
            print >>sys.stderr, '[--- job: prinseq%02d ---] Slurm script written to "prinseq_job.%02d.sh".' % ((i+1), (i+1))
        else:
            # Submit
            from subprocess import Popen,PIPE
            p = Popen(['sbatch','-J','prinseq%02d' % (i+1)],  stdin=PIPE, stdout=PIPE)
            out,err = p.communicate(input='\n'.join(job_sh))
            print >>sys.stderr, '[--- job: prinseq%02d ---] %s. Running %d prinseq processes.' % ((i+1), out.strip('\n'), len(node))

        # If args.detail, print to stderr
        if args.detail:
            print >>sys.stderr, '\n'.join(job_sh)
            print >>sys.stderr, ''

if __name__ == '__main__':
    import argparse
    parser = findfiles.defaultparser('''Filter and quality control for sequencing reads using prinseq''')    
    
    preset_group = parser.add_argument_group("Presets")
    preset_group.add_argument('--illumina', action='store_true',
                             help='''Illumina filtering presets. 
                                     Equivalent to "%s"''' % ' '.join('-%s %s' % (k,v) for k,v in ILLUMINA_PRESET.iteritems()) )
    preset_group.add_argument('--iontorrent', action='store_true',
                             help='''IonTorrent filtering presets. 
                                     Equivalent to "%s"''' % ' '.join('-%s %s' % (k,v) for k,v in IONTORRENT_PRESET.iteritems()) )
    slurm_group = parser.add_argument_group("Slurm options")
    slurm_group.add_argument('--nosubmit', action='store_true',
                             help="Do not submit jobs directly. Output is written to batch scripts")
    slurm_group.add_argument('--walltime', type=int, default=480, 
                             help="Slurm walltime request (in minutes)")
    slurm_group.add_argument('--partition', default="short", 
                             help="Slurm partition request for resource allocation")
    slurm_group.add_argument('--njobs', type=int,
                             help='''Total number of SLURM jobs to start. Each job will use
                                     one node.''')
    slurm_group.add_argument('--nslots', type=int, default=8, 
                             help='''Number of simultaneous processes that can run on a
                                     node. For uncompressed files, this is equal to the 
                                     number of CPUs on each node. For compressed files,
                                     RAM is the limiting factor becausewe are limited by the RAM consumed during decompression
                                     since gzi''')
    slurm_group.add_argument('--nsequential', type=int, default=1, 
                             help='''Number of process to run sequentially on the same node.
                                     If queue wait times are long, it may be advantageous to
                                     run sequentially to avoid waiting for additional
                                     allocations. For example, if you have 18 read sets to
                                     process and nslots is 16, the wrapper will request 2
                                     allocations. However, setting nsequential=2 will
                                     allow jobs to run on the same node after one of the
                                     initial jobs has completed.''')
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
    output_group.add_argument('--outdir',
                             help='''Output directory for files. Default is the same 
                                     directory as input files''')                                      
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
