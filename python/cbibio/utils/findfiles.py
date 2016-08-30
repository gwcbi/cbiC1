#! /usr/bin/env python
import sys
import argparse
from glob import glob
import re
import os
import gzip

from collections import defaultdict

from guess_encoding import guess_encoding


FILE_EXTENSIONS = { 'FASTQ': set(['fastq','fq']),
                    'FASTA': set(['fasta','fa', 'fas', 'fna']),
                    'SRA': set(['sra']),
                  }

EXT2FTYPE = {}
for ft,fxs in FILE_EXTENSIONS.iteritems():
    for fx in list(fxs):
        EXT2FTYPE['.%s' % fx] = ft

'''
EXT2FTYPE =       { '.fastq':  'FASTQ',
                    '.fq':     'FASTQ',
                    '.fasta':  'FASTA',
                    '.fas':    'FASTA',
                    '.fna':    'FASTA',
                    '.fa':     'FASTA',
                    '.sra':    'SRA',
                   }
'''

def sizeof_fmt(num, suffix=''):
    for unit in ['','K','M','G','T']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "fail"


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
        self._encoding = None
    
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
                self._encoding = guess_encoding(gzip.open(self.fullpath), 100)
            else:
                self._encoding = guess_encoding(self.fullpath, 100)
        return self._encoding        
    
    def __str__(self):
        return '\t'.join([self.fullpath, self.filetype, 'gzip' if self.zipped else '', self.dirname, self.basename])

class ReadFileSet:
    ''' Set of read files that belong together, such as paired-end reads '''
    def __init__(self, r1=None, r2=None, up=None):
        self.R1 = r1
        self.R2 = r2
        self.UP = up
    def paired(self):
        return self.R1 is not None and self.R2 is not None and self.UP is None
    def unpaired(self):
        return self.R1 is None and self.R2 is None and self.UP is not None
    def both(self):
        return self.R1 is not None and self.R2 is not None and self.UP is not None
    def status(self):
        if self.paired(): return 'paired'
        if self.unpaired(): return 'unpaired'
        if self.both(): return 'both'
        return 'incomplete'
    def check_pair(self):
        if self.R1 is not None and self.R2 is not None:
            pn1,mn1 = get_read_pairname(self.R1.peek()[0])
            pn2,mn2 = get_read_pairname(self.R2.peek()[0])
            return pn1 == pn2 and (mn1 == '1' or mn1 is None) and (mn2 == '2' or mn2 is None)
        return True
    def size(self):
        ret = self.R1.fsize if self.R1 is not None else 0
        ret += self.R2.fsize if self.R2 is not None else 0
        ret += self.UP.fsize if self.UP is not None else 0
        return ret
    def encoding(self):
        if self.paired() or self.both():
            assert self.R1.encoding() == self.R2.encoding()
            return self.R1.encoding()
        else:
            return self.UP.encoding()
    
    def __str__(self):
        ret = 'ReadFileSet (%s). Total size: %s\n' % (self.status(), sizeof_fmt(self.size()))
        if self.R1 is not None: ret += '\tR1:\t%s\n' % str(self.R1)
        if self.R2 is not None: ret += '\tR2:\t%s\n' % str(self.R2)
        if self.UP is not None: ret += '\tUP:\t%s\n' % str(self.UP)
        return ret.strip('\n')

def get_read_pairname(rn):
    ''' Return the pairname and mate number for given read name '''
    # Illumina <1.8
    m = re.search('^[@>](?P<pairname>[-:\w]+(#[\dATCG]+)?)\s*/(?P<matenum>[123])', rn)
    if m:
        matenum = '2' if m.group('matenum') == '3' else m.group('matenum')
        return m.group('pairname'), matenum
    # Illumina >=1.8         
    m = re.search('^[@>](?P<pairname>[-:\w]+)\s+(?P<matenum>[123]):[:\w]+', rn)
    if m:
        matenum = '2' if m.group('matenum') == '3' else m.group('matenum')
        return m.group('pairname'), matenum
    return rn.split()[0][1:], None

'''
assert get_read_pairname('@HWUSI-EAS100R:6:73:941:1973#0/1') == ('HWUSI-EAS100R:6:73:941:1973#0', '1')
assert get_read_pairname('@HWUSI-EAS100R:6:73:941:1973/1') == ('HWUSI-EAS100R:6:73:941:1973','1')
assert get_read_pairname('@EAS51_105_FC20G7EAAXX_R1:1:1:471:409#ATCACG/2') == ('EAS51_105_FC20G7EAAXX_R1:1:1:471:409#ATCACG','2')
assert get_read_pairname('@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG') == ('EAS139:136:FC706VJ:2:2104:15343:197393', '1')
assert get_read_pairname('@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36') == ('SRR001666.1', None)
assert get_read_pairname('@071112_SLXA-EAS1_s_7:5:1:817:345') == ('071112_SLXA-EAS1_s_7:5:1:817:345', None)
'''

def get_file_pairname(bfn):
    ''' Return the pairname and mate number for given file name '''
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

'''
assert get_file_pairname('P153SEQ0079_P189M0155_S23_L001_R2') == ('P153SEQ0079_P189M0155_S23_L001','2')
assert get_file_pairname('NA10831_ATCACG_L002_R1_001') == ('NA10831_ATCACG_L002','1')
assert get_file_pairname('SRR1551011_pass_1') == ('SRR1551011_pass', '1')
assert get_file_pairname('sample001.1') == ('sample001', '1')
assert get_file_pairname('sample006.R2') == ('sample006', '2')
assert get_file_pairname('SRR1551011_pass') == ('SRR1551011_pass', None)
'''

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

def rec_find_files(start='.'):
    ''' Recursively find files in subfolders of start '''
    for root, dirnames, filenames in os.walk(start):
        for f in filenames:
            yield os.path.join(root,f)

def filter_files(flist, ftype, include_str=None, exclude_str=None):
    ''' Filter list of files according to file type and regular expressions '''
    includeRE = re.compile('') if include_str is None else re.compile(include_str)
    excludeRE = re.compile('\0') if exclude_str is None else re.compile(exclude_str)
    ret = []
    for bf in flist:
        if bf.filetype == ftype and bf.matchRE(includeRE) and bf.notRE(excludeRE):
            ret.append(bf)
    return ret

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
        for rfset in rfsets:
            print rfset
    else:
        rfsets = find_mate_files_byreadname(biofiles)
        for rfset in rfsets:
            print rfset

    print >>sys.stderr, '%d file sets were found:' % (len(rfsets))
    print >>sys.stderr, '\t%d paired sets' % sum(rfset.paired() for rfset in rfsets)
    print >>sys.stderr, '\t%d unpaired sets' % sum(rfset.unpaired() for rfset in rfsets)
    print >>sys.stderr, '\t%d paired+unpaired sets' % sum(rfset.both() for rfset in rfsets)    
