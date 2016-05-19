#! /usr/bin/env python
import sys
from glob import glob

""" Constants """

SCRIPT_HEADER = '''#! /bin/bash
#SBATCH -t %(walltime)s
#SBATCH -p %(partition)s
#SBATCH -N %(nodes)s

SN="fastq-dump-wrapper"
echo "[---$SN---] ($(date)) Starting $SN"
t1=$(date +"%%s")

module load %(sra_module)s
echo "[---$SN---] ($(date))  fastq-dump executable: $(which fastq-dump)"

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

# fastq-dump options:
#     --split-files    Split mate pairs into seperate files
#     --gzip           Compress output using gzip
#     -F               Defline contains only original sequence name
#     -R               Split into files by READ_FILTER value
#     --defline-qual   Defline format specification for quailty. Use "+" to 

FASTQ_DUMP_CMD = '''fastq-dump --split-files --gzip -F -R --defline-qual "+"  %(srafile)s &'''

""" Functions """
def chunks(l, n):
  ''' Yield successive n-sized chunks from list. '''
  for i in xrange(0, len(l), n):
    yield l[i:i+n]

""" main """
def main(args):
    # Get list of files
    if args.fofn is not None:
        allfiles  = sorted([l.strip() for l in args.fofn])
    else:
        allfiles = sorted(glob('*.sra'))
    
    for i,chunk in enumerate(chunks(allfiles, args.chunksize)):
        filecmds = []
        for fs in chunk:
            cmd = FASTQ_DUMP_CMD % {'srafile':fs}
            filecmds.extend(['echo "[---$SN---] ($(date)) COMMAND: %s"' % cmd, cmd, ''])
        
        script = [SCRIPT_HEADER % vars(args)] + filecmds + [SCRIPT_FOOTER % vars(args)]
        if args.nosubmit:
            with open('job.%02d.sh' % (i+1), 'w') as outh: 
                print >>outh, '\n'.join(script)
            print >>sys.stderr, '[--- job: fastq-dump%02d ---] Slurm script written to "job.%02d.sh".' % ((i+1), (i+1))                
        else:
            from subprocess import Popen,PIPE
            p = Popen(['sbatch','-J','fastq-dump%02d' % (i+1)],  stdin=PIPE, stdout=PIPE)
            out,err = p.communicate(input='\n'.join(script))
            # out = 'Dummy slurm response 12345.\n'            
            print >>sys.stderr, '[--- job: fastq-dump%02d ---] %s. Running %d fastq-dump processes.' % ((i+1), out.strip('\n'), len(chunk))            
        
        if args.detail:
            print >>sys.stderr, '\n'.join(script)
            print >>sys.stderr, ''

if __name__ == '__main__':
  import argparse
  import sys
  parser = argparse.ArgumentParser(description='Wrapper for unpacking many SRA files in parallel', 
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  input_group = parser.add_argument_group("Input options")
  input_group.add_argument('--fofn', type=argparse.FileType('r'),
                           help="File of file names")

  slurm_group = parser.add_argument_group("Slurm options")
  slurm_group.add_argument('--nosubmit', action='store_true',
                           help="Do not submit jobs directly. Output is written to batch scripts")
  slurm_group.add_argument('--chunksize', type=int, default=16, 
                           help="Number of instances per job. This should be less than or equal to the number of CPUs per node")
  slurm_group.add_argument('--walltime', type=int, default=480, 
                           help="Slurm walltime request (in minutes)")
  slurm_group.add_argument('--partition', default="short", 
                           help="Slurm partition request for resource allocation")
  slurm_group.add_argument('--nodes', type=int, default=1, 
                           help="Slurm node request")
  slurm_group.add_argument('--sra_module', default="ncbi-sra/2.5.2", 
                           help='Name of sra toolkit module. Calling "module load [sra_module]" must load fastq-dump into environment')
  slurm_group.add_argument('--detail', action='store_true',
                           help='Show detailed information about jobs being submitted.')                           
  args = parser.parse_args()
  main(args)
