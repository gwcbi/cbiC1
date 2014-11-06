#! /usr/bin/env python

class Fastq:
  ''' Fastq record '''
  # Regular expressions for pair suffix
  import re
  re_f = re.compile(r"(\S+)(/1|\.f|\.[sfp]\d\w*)$")
  re_r = re.compile(r"(\S+)(/2|\.r|\.[rq]\d\w*)$")
  re_pair = re.compile(r"(\S+)(/[12]|\.[fr]|\.[sfprq]\d\w*)$")  
  
  def __init__(self,name,seq,plus,qual):
    self.name = name
    self.seq  = seq
    self.plus = plus
    self.qual = qual
  
  def basename(self):
    ''' Returns the read name without pair suffix
        Code here is modified from bitbucket repo peterjc/galaxy-central
        Script is tools/fastq/fastq_paired_unpaired.py
        See https://bitbucket.org/peterjc/galaxy-central

        Cope with three widely used suffix naming convensions,
          Illumina: /1 or /2
          Forward/revered: .f or .r
          Sanger, e.g. .p1k and .q1k
          See http://staden.sourceforge.net/manual/pregap4_unix_50.html
    '''
    m = self.re_pair.search(self.name)
    if m: return m.group(1)
    return self.name    

  def orient(self):
    if self.re_f.search(self.name): return 'fwd'
    if self.re_r.search(self.name): return 'rev'
    return 'no'
    
  def __str__(self):
    return '\n'.join([self.name,self.seq,self.plus,self.qual])

class FastqFile:
  def __init__(self,f):
    self.fh = open(f,'rU') if isinstance(f,basestring) else f
  
  def records(self):
    while self.fh:
      yield Fastq(*[self.fh.next().strip('\n') for _ in range(4)])

"""
from cbibio.utils.Fastq import Fastq
template_name = '@B206TFABXX:6:67:12897:87553'
suffixes = ['.f','.r','/1','/2','.p1','.p1k','.p1lk','.q1','.q1k','.q1lk','.s1','.s1k','']
for suf in suffixes:
  print Fastq('%s%s' % (template_name,suf),'ATG','+','###').basename()

for suf in suffixes:
  print '%s%s -> %s' % (template_name,suf,Fastq('%s%s' % (template_name,suf),'ATG','+','###').orient())
  
  assert Fastq('%s%s' % (template_name,suf),'ATG','+','###').basename() == template_name
"""
