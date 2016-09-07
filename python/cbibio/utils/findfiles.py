#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Functions and classes for working with the filesystem


Attributes:
    FILE_EXTENSIONS (:obj:`dict`): Dictionary mapping file types to a set of expected (or 
        acceptable) extensions
    EXT2FTYPE (:obj:`dict`): Dictionary mapping file extensions to the file type
"""
__author__ = "Matthew L. Bendall"

import sys
import argparse
from glob import glob
import re
import os
import gzip

from collections import defaultdict
from cbibio.utils.guess_encoding import guess_encoding

""" Constants """

FILE_EXTENSIONS = { 'FASTQ': set(['fastq','fq']),
                    'FASTA': set(['fasta','fa', 'fas', 'fna']),
                    'SRA': set(['sra']),
                  }

EXT2FTYPE = {}
for ft,fxs in FILE_EXTENSIONS.iteritems():
    for fx in list(fxs):
        EXT2FTYPE['.%s' % fx] = ft


""" Functions working with files """

def fmtsize_h(num, suffix=''):
    """ Format file size so it is human readable

    Args:
        num (int): File size in bytes
        suffix (str, optional): Add suffix to end of string

    Returns:
        str: Human readable file size        
    """
    for unit in ['','K','M','G','T']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "fail"

def recursively_find_files(top='.'):
    """ Recursively generate files in directory tree

    Args:
        top (str): The root of the directory tree (inclusive)
        
    Yields:
        str: The next file in the directory tree
    """
    for root, dirnames, filenames in os.walk(top=top):
        for f in filenames:
            fullpath = os.path.join(root, f)
            if os.path.isfile(fullpath):
                yield fullpath

def filter_files(flist, ftype=None, include=None, exclude=None):
    """ Filter files by filename and filetype
    
    If both include and exclude are used, file must meet both criteria (and, not or)
    
    Args:
        flist (:obj:`list` of :obj:`BioFile`): List of BioFile to filter
        ftype (str or list of str, optional): File type(s) to include
        include (str, optional): File will be included only if filename matches regex
        exclude (str, optional): File will be excluded if filename matches regex

    Returns:
        :obj:`list` of :obj:`BioFile`: List of files that meet filtering criteria
    """
    ret = []
    if type(ftype) is str: ftype = [ ftype ]
    for bf in flist:
        cond = True if ftype is None else bf.filetype in ftype
        cond &= True if include is None else bf.matchRE(re.compile(include))
        cond &= True if exclude is None else bf.notRE(re.compile(exclude))
        if cond:
            ret.append(bf)
    return ret

def get_read_pairname(rn):
    """ Get the readname and mate number for given readname
    
    Examples:
        >>> get_read_pairname('@HWUSI-EAS100R:6:73:941:1973#0/1')
        ('HWUSI-EAS100R:6:73:941:1973#0', '1')
        >>> get_read_pairname('@HWUSI-EAS100R:6:73:941:1973/1')
        ('HWUSI-EAS100R:6:73:941:1973','1')
        >>> get_read_pairname('@EAS51_105_FC20G7EAAXX_R1:1:1:471:409#ATCACG/2')
        ('EAS51_105_FC20G7EAAXX_R1:1:1:471:409#ATCACG','2')
        >>> get_read_pairname('@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG')
        ('EAS139:136:FC706VJ:2:2104:15343:197393', '1')
        >>> get_read_pairname('@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36')
        ('SRR001666.1', None)
        >>> get_read_pairname('@071112_SLXA-EAS1_s_7:5:1:817:345')
        ('071112_SLXA-EAS1_s_7:5:1:817:345', None)
        >>> get_read_pairname('@FC61B8BAAXX:4:1:1:174#TNACCA/1')
        ('@FC61B8BAAXX:4:1:1:174#TNACCA', '1')
    """
    # Illumina <1.8
    m = re.search('^[@>](?P<pairname>[-:\w]+(#[\dA-Z]+)?)\s*/(?P<matenum>[123])', rn)
    if m:
        matenum = '2' if m.group('matenum') == '3' else m.group('matenum')
        return m.group('pairname'), matenum
    # Illumina >=1.8
    m = re.search('^[@>](?P<pairname>[-:\w]+)\s+(?P<matenum>[123]):[:\w]+', rn)
    if m:
        matenum = '2' if m.group('matenum') == '3' else m.group('matenum')
        return m.group('pairname'), matenum
    return rn.split()[0][1:], None

def get_file_pairname(bfn):
    """ Return the pairname and mate number for given file name
    
    Examples:
        >>> get_file_pairname('P153SEQ0079_P189M0155_S23_L001_R2')
        ('P153SEQ0079_P189M0155_S23_L001','2')
        >>> get_file_pairname('NA10831_ATCACG_L002_R1_001')
        ('NA10831_ATCACG_L002','1')
        >>> get_file_pairname('SRR1551011_pass_1')
        ('SRR1551011_pass', '1')
        >>> get_file_pairname('sample001.1')
        ('sample001', '1')
        >>> get_file_pairname('sample006.R2')
        ('sample006', '2')
        >>> get_file_pairname('SRR1551011_pass')
        ('SRR1551011_pass', None)
    """
    # Example: 'P153SEQ0079_P189M0155_S23_L001_R2'
    m = re.search('(?P<pairname>[\w\.]+)[\._]R?(?P<matenum>[12])$', bfn)
    if m:
        return m.group('pairname'), m.group('matenum')
    # Illumina FASTQ files use the following naming scheme:
    # <sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz        
    # Example: bfn = 'NA10831_ATCACG_L002_R1_001'
    m = re.search('(?P<pairname>[\w\.]+)_R(?P<matenum>[12])_\d{3}', bfn)
    if m:
        return m.group('pairname'), m.group('matenum')
    return bfn, None

def find_mate_files_byfilename(flist, check_readnames=True):
    ''' Identify paired files by examining the filenames '''
    # First pass: find files with matching pairnames
    pass1 = defaultdict(list)
    for bf in flist:
        pn,mn = get_file_pairname(bf.basename)
        pass1[pn].append((bf,mn))
    
    # Second pass: check number of files and mate numbers
    pass2 = []
    for pn,tups in pass1.iteritems():
        if len(tups)==1:
            # Unpaired sample
            pass2.append(ReadFileSet(up=tups[0][0]))
        elif len(tups)==2:
            if set(t[1] for t in tups) == set(['1', '2']):
                r1,r2 = [t[0] for t in sorted(tups,key=lambda x:int(x[1]))]
                pass2.append(ReadFileSet(r1=r1,r2=r2))
            else:
                print >>sys.stderr, 'Pairing Error:\n%s' % '\n'.join(str(t[0]) for t in tups)
                pass2.extend([ReadFileSet(up=t[0]) for t in tups])
        elif len(tups)==3:
            if set(t[1] for t in tups) == set(['1', '2', None]):
                tmp = {t[1]:t[0] for t in tups}
                pass2.append(ReadFileSet(r1=tmp['1'], r2=tmp['2'], up=tmp[None]))
            else:
                print >>sys.stderr, 'Pairing Error:\n%s' % '\n'.join(str(t[0]) for t in tups)
                pass2.extend([ReadFileSet(up=t[0]) for t in tups])
        else:
            print >>sys.stderr, 'Pairing Error:\n%s' % '\n'.join(str(t[0]) for t in tups)
            pass2.extend([ReadFileSet(up=t[0]) for t in tups])
    
    if not check_readnames:
        print >>sys.stderr, 'Bypassing check'
        return pass2
    
    # Final pass: check that read names match
    checked_sets = []
    for rfset in pass2:
        if rfset.check_pair():
            checked_sets.append(rfset)
        else:
            print >>sys.stderr, 'PAIRING ERROR: Files seem to match, but read names do not match.'
            print >>sys.stderr, '               File1: %s\tReadname1: %s' % (rfset.R1.fullpath, rfset.R1.peek()[0])
            print >>sys.stderr, '               File2: %s\tReadname2: %s' % (rfset.R2.fullpath, rfset.R2.peek()[0])
            checked_sets.append(ReadFileSet(up=rfset.R1))
            checked_sets.append(ReadFileSet(up=rfset.R2))
            if rfset.UP is not None:
                checked_sets.append(ReadFileSet(up=rfset.UP))
    
    return checked_sets

def find_mate_files_byreadname(flist):
    ''' Identify paired files by examining the first readname '''
    pairs = defaultdict(list)
    for bf in flist:
        pn,mn = get_read_pairname(bf.peek()[0])
        pairs[pn].append((bf,mn))
    
    ret = []
    for pn,tups in pairs.iteritems():
        if len(tups)==1:
            # Unpaired sample
            ret.append(ReadFileSet(up=tups[0][0]))
        elif len(tups)==2:
            if set(t[1] for t in tups) == set(['1', '2']):
                r1,r2 = [t[0] for t in sorted(tups,key=lambda x:int(x[1]))]
                ret.append(ReadFileSet(r1=r1,r2=r2))
            else:
                # The mate numbers are incorrect or not included in read name.
                # Instead determine R1 and R2 by file basename
                r1,r2 = [t[0] for t in sorted(tups,key=lambda x:x[0].basename)]
                ret.append(ReadFileSet(r1=r1,r2=r2))                
        elif len(tups)==3:
            if set(t[1] for t in tups) == set(['1', '2', None]):
                tmp = {t[1]:t[0] for t in tups}
                ret.append(ReadFileSet(r1=tmp['1'], r2=tmp['2'], up=tmp[None]))
            else:
                print >>sys.stderr, 'Pairing Error:\n%s' % '\n'.join(str(t[0]) for t in tups)
                ret.extend([ReadFileSet(up=t[0]) for t in tups])
        else:
            print >>sys.stderr, 'Pairing Error:\n%s' % '\n'.join(str(t[0]) for t in tups)
            ret.extend([ReadFileSet(up=t[0]) for t in tups])
    
    return ret

def find_mate_files(flist, method, check_readnames):
    if method == 'readname':
        return find_mate_files_byreadname(flist)
    if method == 'filename':
        return find_mate_files_byfilename(flist, check_readnames)
    assert False, "ERROR: unknown method for finding mate files - %s" % method

class BioFile:
    ''' Store metadata about a file '''
    def __init__(self, fpath):
        self.fullpath = fpath    
        root,ext = os.path.splitext(fpath)
        self.zipped = ext.lower() == '.gz'
        if self.zipped:
            root,ext = os.path.splitext(root)
        self.dirname, self.basename = os.path.split(root)
        self.filetype = EXT2FTYPE[ext.lower()] if ext.lower() in EXT2FTYPE else 'unk'
        self.fsize = os.path.getsize(self.fullpath)
        self._filecache = []
        self._encoding = None # Determined the first time self.encoding() is called
    
    def matchRE(self, regex):
        # Return True if basename matches regex
        return regex.search(self.basename) is not None
    
    def notRE(self, regex):
        # Return True if basename does not match regex    
        return regex.search(self.basename) is None
    
    def peek(self, nlines=1):
        if len(self._filecache) >= nlines:
            return self._filecache[:nlines]
        if len(self._filecache)>0:
            # We've done this before
            print >>sys.stderr, "Going to disk"
        if self.zipped:
            lines = (l.strip('\n') for l in gzip.open(self.fullpath, 'rb'))
        else:
            lines = (l.strip('\n') for l in open(self.fullpath, 'r'))
        self._filecache = [lines.next() for _ in range(nlines)]
        return self._filecache #[l for l in range(nlines)]
    
    def encoding(self):
        if self.filetype is 'FASTQ' and self._encoding is None:
            if self.zipped:
                self._encoding = guess_encoding(gzip.open(self.fullpath, 'rb'), 100)
            else:
                self._encoding = guess_encoding(self.fullpath, 100)
        return self._encoding
    
    def __str__(self):
        return '\t'.join([self.fullpath, self.filetype, 'gzip' if self.zipped else '', self.dirname, self.basename])

class ReadFileSetError(Exception):
    """ Raised when ReadFileSet does not meet expectations """
    pass

class ReadFileSet:
    ''' Set of read files that belong together, such as paired-end reads '''
    def __init__(self, r1=None, r2=None, up=None):
        self.R1 = r1
        self.R2 = r2
        self.UP = up
        self._set_status()
        self._encoding = None # Determined the first time self.encoding() is called
    
    def _set_status(self):
        """ Determine attributes for ReadFileSet """
        self.paired   = (self.R1 is not None) & (self.R2 is not None) & (self.UP is None)
        if self.paired:
            if self.R1.filetype != self.R2.filetype:
                msg = "ERROR (ReadFileSetError): Files in pair must have same format"
                raise ReadFileSetError('%s\n  %s\n  %s' % (msg, self.R1,self.R2))            
            if self.R1.zipped != self.R2.zipped:
                msg = "ERROR (ReadFileSetError): Files in pair must have same compression"
                raise ReadFileSetError('%s\n  %s\n  %s' % (msg, self.R1,self.R2))
            self.filetype = self.R1.filetype
            self.zipped = self.R1.zipped
            self.status = 'paired'
        
        self.unpaired = (self.R1 is None) & (self.R2 is None) & (self.UP is not None)
        if self.unpaired:
            self.filetype = self.UP.filetype
            self.zipped = self.UP.zipped
            self.status = 'unpaired'
        
        self.mixed    = (self.R1 is not None) & (self.R2 is not None) & (self.UP is not None)
        if self.mixed:
            if self.R1.filetype != self.R2.filetype or self.R1.filetype != self.UP.filetype:
                msg = "ERROR (ReadFileSetError): Files in set must have same format"
                raise ReadFileSetError('%s\n  %s\n  %s\n  %s' % (msg, self.R1, self.R2, self.UP))
            if self.R1.zipped != self.R2.zipped or self.R1.zipped != self.UP.zipped:
                msg = "ERROR (ReadFileSetError): Files in set must have same compression"
                raise ReadFileSetError('%s\n  %s\n  %s\n  %s' % (msg, self.R1, self.R2, self.UP))
            self.filetype = self.R1.filetype
            self.zipped = self.R1.zipped
            self.status = 'mixed'
    
    def check_pair(self):
        """ Check that pairing is correct
        
            This function compares the first read name of each file to see whether they
            match.
        """
        if self.R1 is None or self.R2 is None:
            return True # Not paired, so return True
        pn1,mn1 = get_read_pairname(self.R1.peek()[0])
        pn2,mn2 = get_read_pairname(self.R2.peek()[0])
        return pn1 == pn2 and (mn1 == '1' or mn1 is None) and (mn2 == '2' or mn2 is None)
    
    def size(self, h=False):
        """ Returns sum of file sizes for set
        Args:
            h (bool): Return a human readable string
        
        Returns:
            int: Sum of all file sizes. If h is True, returns human-readable string.
        """
        if self.paired:
            ret = self.R1.fsize + self.R2.fsize
        elif self.mixed:
            ret = self.R1.fsize + self.R2.fsize + self.UP.fsize
        elif self.unpaired:
            ret = self.UP.fsize
        return fmtsize_h(ret) if h else ret
        
    def encoding(self):
        """ Determine the encoding for FASTQ files
        
        Returns:
            str: Either 'Phred+33' or 'Phred+64' depending on encoding offset. Defaults to
                None if files are not fastq
        """
        if self.filetype is 'FASTQ' and self._encoding is None: # Encoding has not been determined
            if self.paired:
                if self.R1.encoding() != self.R2.encoding():
                    msg = "ERROR (ReadFileSetError): Files in pair must have same encoding"
                    raise ReadFileSetError('%s\n  %s\n  %s' % (msg, self.R1, self.R2))
                self._encoding = self.R1.encoding()
            elif self.mixed:
                if self.R1.encoding() != self.R2.encoding() or self.R1.encoding() != self.UP.encoding():
                    msg = "ERROR (ReadFileSetError): Files in set must have same encoding"
                    raise ReadFileSetError('%s\n  %s\n  %s\n  %s' % (msg, self.R1, self.R2, self.UP))
                self._encoding = self.R1.encoding()
            elif self.unpaired:
                self._encoding = self.UP.encoding()
        return self._encoding

    def output_prefix(self):
        if self.unpaired:
            return self.UP.dirname, self.UP.basename
        elif self.paired or self.mixed:
            pn1,mn1 = get_file_pairname(self.R1.basename)
            pn2,mn2 = get_file_pairname(self.R2.basename)
            if pn1 == pn2 and self.R1.dirname == self.R2.dirname:
                return self.R1.dirname, pn1
            else:
                print >>sys.stderr, "WARNING: Different output prefix:\n  %s\n  %s" % (self.R1, self.R2)
                return self.R1.dirname, pn1
    
    def __str__(self):
        ret = 'ReadFileSet (%s). Total size: %s\n' % (self.status, self.size(h=True))
        if self.R1 is not None: ret += '\tR1:\t%s\n' % str(self.R1)
        if self.R2 is not None: ret += '\tR2:\t%s\n' % str(self.R2)
        if self.UP is not None: ret += '\tUP:\t%s\n' % str(self.UP)
        return ret.strip('\n')


def find_readsets(infiles, ftypes=['FASTQ','FASTA'], include=None, exclude=None,
                  pair_method="filename", check_readnames=True, single=False, verbose=True):
    if type(infiles) is str:
        biofiles = [BioFile(f) for f in recursively_find_files(infiles)]
    else:
        biofiles = [BioFile(f) for f in infiles if os.path.isfile(f)]
    
    if len(biofiles) == 0:
        if verbose: print >>sys.stderr, "No files found"
        return []
    
    # Filter the file list
    biofiles = filter_files(biofiles, ftypes, include, exclude)
    
    if len(biofiles) == 0:
        if verbose: print >>sys.stderr, "No files remain after filter"
        return []
    
    if verbose:
        print >>sys.stderr, '%d read files were found' % (len(biofiles))
        for bf in biofiles:
            print >>sys.stderr, bf

    if single:
        print >>sys.stderr, "Assuming all files are unpaired"
        rfsets = [ReadFileSet(up=bf) for bf in biofiles]
    else:       
        if verbose:
            s = "Finding mate files using %s" % pair_method
            if pair_method == 'filename':
                if check_readnames: s += ' and checking readnames'
                else: s += ' without checking readnames'
            print >>sys.stderr, s
    
        rfsets = find_mate_files(biofiles, pair_method, check_readnames)
    
    if verbose:
        for rfset in rfsets:
            print >>sys.stderr, rfset
        print >>sys.stderr, '%d file sets were found:' % (len(rfsets))
        print >>sys.stderr, '\t%d paired sets' % sum(rfset.paired for rfset in rfsets)
        print >>sys.stderr, '\t%d unpaired sets' % sum(rfset.unpaired for rfset in rfsets)
        print >>sys.stderr, '\t%d paired+unpaired sets' % sum(rfset.mixed for rfset in rfsets)

    return rfsets

def defaultparser(desc=''):
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    input_group = parser.add_argument_group("Input options")
    input_group.add_argument('--include_files',
                        help='''Only include files that match the given search string. Can
                                be a regular expression.''')
    input_group.add_argument('--exclude_files',
                        help='''Exclude files that match the given search string. Can be a
                                regular expression.''')
    input_group.add_argument('--file_type', default="FASTQ", choices=FILE_EXTENSIONS.keys(),
                        help='''Types of files to find.''')
    input_group.add_argument('--single', action='store_true',
                        help='''Assume all files are single end. Do not attempt to find pairs''')
    input_group.add_argument('--by_readname', action='store_true',
                        help='''Use readname instead of filename to find pairs. The
                                advantage here is that the files can have any random name,
                                as long as the read names match they will be paired. There
                                are some pitfalls that might cause this to fail: 1) If you
                                have multiple copies of the same file with different
                                filenames; and 2) If some file sets include both paired and 
                                unpaired files.''')
    input_group.add_argument('--no_check', action='store_true',
                        help='''Do not check the readnames when finding pairs by filename.''')
    input_group.add_argument('--quiet', action='store_true',
                        help='''Reporting mode''')
    input_group.add_argument('--fofn', type=argparse.FileType('r'),
                        help='''Provide a file of file names instead of infiles''')
    input_group.add_argument('infiles', nargs='*', default=".",
                             help='''Input files. Can be a list of paths for each file or
                                     a wildcard expression (handled by the OS). 
                                     Recursively searches current directory if no argument
                                     is present. ''')
    return parser

def promptabort():
    response = None
    while response is None:
        response = raw_input("Do you wish to continue? [Y/N] ")
        if response in set(['N','n','no','NO','No']):
            sys.exit("Aborted")
        elif response not in set(['Y','y','yes','Yes','YES']):
            response = None
    assert response in set(['Y','y','yes','Yes','YES'])

if __name__=='__main__':
    parser = defaultparser('Run a test of file discovery scripts.')
    args = parser.parse_args()
    if args.fofn:
        infiles = [l.strip('\n') for l in args.fofn if l.strip()]
    else:
        infiles = args.infiles

    pair_method = 'readname' if args.by_readname else 'filename'
    find_readsets(infiles, args.file_type, args.include_files, args.exclude_files,
                  pair_method, not args.no_check, args.single, not args.quiet)

