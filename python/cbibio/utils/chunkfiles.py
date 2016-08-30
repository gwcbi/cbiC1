#! /usr/bin/env python

from findfiles import *
import math
import sys

def distribute_readsets(rfsets, nslots, nsequential=1, nnodes=None):
    '''
    Args:
        nslots(int): number of simultaneous processes that can be run on a node. This is
                     provided as -j to parallel
        nnodes(int): total number of nodes being requested
    '''
    nproc = len(rfsets) # The number of processes that need to run
    if nnodes is None:
        max_proc_per_node = nslots * nsequential # Maximum number of processes that will run on 1 node
        nnodes = (nproc / max_proc_per_node) + (nproc % max_proc_per_node > 0)
    else:
        if nnodes * nslots * nsequential < nproc:
            nsequential = (nproc / (nnodes * nslots)) + ((nproc % (nnodes * nslots)) > 0)
        max_proc_per_node = nslots * nsequential
    
    rfsets.sort(key=lambda x:-x.size())
    nodes = [list() for _ in range(nnodes)]
    for j,rfset in enumerate(rfsets):
        # Option 1: Assign process to the smallest node, until maximum
        # node_sizes = sorted([(i,sum(_.size() for _ in n)) for i,n in enumerate(nodes)], key=lambda x:x[1])
        # destnode = [i for i,sz in node_sizes if len(nodes[i]) < max_proc_per_node][0]
        # Option 2: Assign process sequentially, i.e. 0,1,2,3,0,1,2,3,0,1,2. 
        # Node 0 will end up with the most data
        # destnode = j%len(nodes)
        # Option 3: Assign process ascending/descending, i.e. 0,1,2,3,3,2,1,0,0,1,2
        # More balanced than option 2, doesn't require size like option 1
        nn = len(nodes)
        destnode = (j%nn) if (j/nn)%2==0 else (nn-(j%nn)-1)
        # Add process to destination node
        nodes[destnode].append(rfset)
    
    return nodes

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Run a test of file discovery scripts.')
    parser.add_argument('--include_files',
                        help='''Only include files that match the given search string. Can
                                be a regular expression.''')
    parser.add_argument('--exclude_files',
                        help='''Exclude files that match the given search string. Can be a
                                regular expression.''')
    parser.add_argument('--file_type', default="FASTQ", choices=FILE_EXTENSIONS.keys(),
                        help='''Types of files to find.''')
    parser.add_argument('--by_readname', action='store_true',
                        help='''Use readname instead of filename to find pairs. The
                                advantage here is that the files can have any random name,
                                as long as the read names match they will be paired. There
                                are some pitfalls that might cause this to fail: 1) If you
                                have multiple copies of the same file with different
                                filenames; and 2) If some file sets include both paired and 
                                unpaired files.''')
    parser.add_argument('--no_check', action='store_true',
                        help='''Do not check the readnames when finding pairs by filename.''')
    parser.add_argument('--nslots', type=int, default=16,
                        help='''Number of simultaneous processes that can be run on a
                                node. For single-threaded programs, this is equal to the
                                number of CPUs on each node. Otherwise, nslots should be
                                the total number of CPUs divided by the number of CPUs
                                used by each process (rounded up). For example, if each
                                process uses 300% CPU and each node has 16 CPUs, nslots
                                should be set to 5''')
    parser.add_argument('--nnodes', type=int,
                        help='''Number of nodes to run on.''')
    parser.add_argument('--nsequential', type=int, default=1,
                        help='''Number of processes to run sequentially.''')                        
    parser.add_argument('infiles', nargs='*', default=".")
    
    args = parser.parse_args()
    if args.infiles == '.':
        biofiles = [BioFile(f) for f in rec_find_files('.') if os.path.isfile(f)]
    else:
        biofiles = [BioFile(f) for f in args.infiles if os.path.isfile(f)]

    if len(biofiles) == 0: sys.exit("No matching files were found.")
    biofiles = filter_files(biofiles, args.file_type, args.include_files, args.exclude_files)
    if len(biofiles) == 0: sys.exit("No matching files were found.")
    
    print >>sys.stderr, '%d %s files were found' % (len(biofiles), args.file_type)

    if args.file_type != 'FASTQ' and args.file_type != 'FASTA':
        for bf in biofiles:
            print >>sys.stderr, bf
        sys.exit()
    
    print >>sys.stderr, 'Finding mate files'
    if not args.by_readname:
        rfsets = find_mate_files_byfilename(biofiles, not args.no_check)
    else:
        rfsets = find_mate_files_byreadname(biofiles)
    
    print >>sys.stderr, '%d file sets were found:' % (len(rfsets))
    print >>sys.stderr, '\t%d paired sets' % sum(rfset.paired() for rfset in rfsets)
    print >>sys.stderr, '\t%d unpaired sets' % sum(rfset.unpaired() for rfset in rfsets)
    print >>sys.stderr, '\t%d paired+unpaired sets' % sum(rfset.both() for rfset in rfsets)    
    
    nodes =  distribute_readsets(rfsets, args.nslots, args.nsequential, args.nnodes)
    for i,node in enumerate(nodes):
        print >>sys.stderr, '#################################### Node%02d ####################################' % i
        s = '# nproc: %d size: %s' % (len(node), sizeof_fmt(sum(_.size() for _ in node)))
        print >>sys.stderr, '%s#' % s.ljust(79)
        print >>sys.stderr, '################################################################################' 
        print '\n'.join(str(rfset) for rfset in node)
