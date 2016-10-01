#! /usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'Matthew L. Bendall'

import sys
import os

from cbibio.utils import findfiles
from cbibio.utils import chunkfiles

""" Constants """

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
#SBATCH -N 1

SN="bowtie2-wrapper"
echo "[---$SN---] ($(date)) Starting $SN"
t1=$(date +"%%s")

module load cbiC1
module load %(bowtie2_module)s
module load %(samtools_module)s
echo "[---$SN---] ($(date))  bowtie2 executable: $(which bowtie2)"
echo "[---$SN---] ($(date))  samtools executable: $(which samtools)"

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

def fileset_command(rfs, args):
    pipeline = []
    # Pipe in sequence data
    if rfs.paired:
        pipeline.append('seqtk mergepe %s %s' % (rfs.R1.fullpath, rfs.R2.fullpath))
    else:
        if rfs.zipped:
            pipeline.append('gunzip -c %s' % rfs.UP.fullpath)
        else:
            pipeline.append('cat %s' % rfs.UP.fullpath)
    
    # Adjust encoding if necessary (seqtk)
    if rfs.encoding() == 'Phred+64':
        pipeline.append('seqtk seq -Q64 -V -')
    
    # End trimming (seqtk)
    if args.trim_left or args.trim_right:
        nleft  = args.trim_left if args.trim_left else 0
        nright = args.trim_right if args.trim_right else 0        
        pipeline.append('seqtk trimfq -b %d -e %d -' % (nleft, nright))
    
    # Quality trimming (seqtk)
    if args.trim_qual:
        pipeline.append('seqtk trimfq -q %f -' % (args.trim_qual))
    
    # Filter by length or ambiguous characters (seqtk)
    argC = '' if args.keep_comment else '-C'
    argL = '-L %d' % args.min_len if args.min_len else ''
    argN = '-N' if args.drop_ambig else ''
    pipeline.append(' '.join(['seqtk', 'seq', argC, argL, argN, '-']))
    
    # Filter by percent of high-quality bases (fastx)
    if args.fastx_qual:
        pipeline.append('fastq_quality_filter -Q33 -q %d -p %d' % (args.fastx_qual, args.fastx_qual_pct))

    # Filter artifacts (fastx)
    if args.fastx_artifacts:
        pipeline.append('fastx_artifacts_filter -Q33')
    
    # Output   
    indir, bn = rfs.output_prefix()
    if args.outdir:
        prefix = os.path.join(args.outdir, bn)
    else:
        prefix = os.path.join(indir, bn)
    
    if rfs.paired:
        out1 = '%s.seqtk_1.fq' % prefix
        out2 = '%s.seqtk_2.fq' % prefix       
        # Drop pairs where one mate is filtered out (seqtk)
        pipeline.append('seqtk dropse')
        # Split to two files
        pipeline.append('deinterleave_fastq %s %s' % (out1,out2))
        return ' | '.join(pipeline)
    else:
        out1 = '%s.seqtk.fq' % prefix
        return '%s > %s' % (' | '.join(pipeline), out1)

def main(args):
    raise NotImplementedError("Implementation for Bowtie2 wrapper is in progress")
    
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

        body = []
        for rfs in node:
            body.append(fileset_command(rfs, args))

        body.append('EOF')
        job_sh = [header] + body + [footer]
        
        if args.nosubmit:
            with open('job.%02d.sh' % (i+1), 'w') as outh: 
                print >>outh, '\n'.join(job_sh)
            print >>sys.stderr, '[--- job: seqtk%02d ---] Slurm script written to "job.%02d.sh".' % ((i+1), (i+1))
        else:
            # Submit
            from subprocess import Popen,PIPE
            p = Popen(['sbatch','-J','seqtk%02d' % (i+1)],  stdin=PIPE, stdout=PIPE)
            out,err = p.communicate(input='\n'.join(job_sh))
            print >>sys.stderr, '[--- job: seqtk%02d ---] %s. Running %d seqtk processes.' % ((i+1), out.strip('\n'), len(node))
        
        # If args.detail, print to stderr
        if args.detail:
            print >>sys.stderr, '\n'.join(job_sh)
            print >>sys.stderr, ''

if __name__ == '__main__':
    import argparse
    parser = findfiles.defaultparser('''Filter and quality control for sequencing reads using seqtk and fastx toolkit''')
    
    preset_group = parser.add_argument_group("Presets")
    preset_group.add_argument('--telescope', action='store_true',
                              help='''Alignment settings optimized for telescope.
                                    Equivalent to "%s"''' % ' '.join('--%s %s' % (k,v) for k,v in ILLUMINA_PRESET.iteritems()) )
    preset_group.add_argument('--pathoscope', action='store_true',
                              help='''Alignment settings optimized for telescope.
                                    Equivalent to "%s"''' % ' '.join('--%s %s' % (k,v) for k,v in ILLUMINA_PRESET.iteritems()) )
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
    slurm_group.add_argument('--nslots', type=int, default=1, 
                             help='''Number of simultaneous processes that can run on a
                                     node. Since we are running bowtie2 
                                     baseline of 2 CPUs. Running fastx_qual or fastx_artifacts
                                     each require an additional CPU. This means that 4 
                                     processes including both flags can run simultaneously
                                     on a 16-core node.''')
    slurm_group.add_argument('--nsequential', type=int, default=1, 
                             help='''Number of process to run sequentially on the same node.
                                     If queue wait times are long, it may be advantageous to
                                     run sequentially to avoid waiting for additional
                                     allocations. For example, if you have 5 read sets to
                                     process and nslots is 4, the wrapper will request 2
                                     allocations. However, setting nsequential=2 will
                                     allow jobs to run on the same node after an initial
                                     job has completed.''')
    slurm_group.add_argument('--bowtie2_module', default="bowtie2/2.2.3", 
                             help='Name of bowtie2 module. Calling "module load [bowtie2_module]" must load bowtie2 into environment')
    slurm_group.add_argument('--samtools_module', default="samtools/1.3", 
                             help='Name of samtools module. Calling "module load [samtools_module]" must load samtools into environment')
    slurm_group.add_argument('--detail', action='store_true',
                             help='Show detailed information about jobs being submitted.')
    
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument('--keep_comment', action='store_true',
                             help='''Keep comments at the header lines''')
    output_group.add_argument('--outdir',
                             help='''Output directory for files. Default is the same 
                                     directory as input files''')
    
    filter_group = parser.add_argument_group("Filter options")
    filter_group.add_argument('--min_len', type=int,
                              help=" Filter sequence shorter than min_len.")
    filter_group.add_argument('--drop_ambig', action='store_true',
                              help=" Drop sequences containing ambiguous bases.")
    filter_group.add_argument('--min_qual_score', 
                              help='''Filter sequence with at least one quality score below min_qual_score.''')                                  
    
    trim_group = parser.add_argument_group("Trim options")
    trim_group.add_argument('--trim_left', type=int,
                            help="Trim sequence at the 5'-end by trim_left positions.")
    trim_group.add_argument('--trim_right', type=int,
                            help="Trim sequence at the 3'-end by trim_right positions.")
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
