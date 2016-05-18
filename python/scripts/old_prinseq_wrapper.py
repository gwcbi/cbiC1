#! /usr/bin/env python

import argparse
import sys
from subprocess import Popen, PIPE

p = Popen('prinseq-lite -h', shell=True, stdout=PIPE, stderr=PIPE)
sout,serr = p.communicate()
prinseq_help = sout.split('\n')

opts = {}
option_groups = ['OVERALL OPTIONS']
opts[option_groups[-1]] = []
opt = None

for l in prinseq_help[prinseq_help.index('Options:')+1:]:
    if l.strip().startswith('***'):
        if opt is not None:
            opts[option_groups[-1]].append((opt,val, helptext.replace('%','pct.')))
            opt = None
        option_groups.append(l.strip().strip('*').strip())
        opts[option_groups[-1]] = []
    elif l.startswith('    -'):
        if opt is not None:
            opts[option_groups[-1]].append((opt,val, helptext.replace('%','pct.')))
            # print '\n%s\t%s\n%s\n' % (opt,val, helptext)
        opt = l.strip().split()[0]
        val = ''
        helptext = ''
        if l.strip().split()[-1].startswith('<'):
            val = l.strip().split()[-1]
    else:
        if not l.strip():
            helptext += ''
        else:
            helptext += ' %s' % l.strip()

# print '\n'.join(str(_) for _ in opts['FILTER OPTIONS'])

# sys.exit()

parser = argparse.ArgumentParser(description='') # , formatter_class=argparse.ArgumentDefaultsHelpFormatter)

for optgroup in option_groups[1:-1]:
    gparser = parser.add_argument_group(optgroup)
    for opt, val, helptext in opts[optgroup]:
        gparser.add_argument('-%s' % opt, help='%s' % helptext)

"""
for opt, val, helptext in opts['OUTPUT OPTIONS']:
        parser.add_argument('-%s' % opt, help='%s' % helptext)
        
for opt, val, helptext in opts['FILTER OPTIONS']:
        parser.add_argument('-%s' % opt, help=str(helptext))

for opt, val, helptext in opts['TRIM OPTIONS']:
        parser.add_argument('-%s' % opt, help='%s' % helptext)

for opt, val, helptext in opts['REFORMAT OPTIONS']:
        parser.add_argument('-%s' % opt, help='%s' % helptext)

for opt, val, helptext in opts['SUMMARY STATISTIC OPTIONS']:
        parser.add_argument('-%s' % opt, help='%s' % helptext)
"""

# OUTPUT OPTIONS
# FILTER OPTIONS
# TRIM OPTIONS
# REFORMAT OPTIONS
# SUMMARY STATISTIC OPTIONS

parser.parse_args()
