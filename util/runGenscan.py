#!/usr/bin/env python

# Copyright (c) 2006
# Colin Dewey (University of Wisconsin-Madison)
# cdewey@biostat.wisc.edu
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


import sys
import os

import FASTA
import GFF

from optparse import OptionParser

usage = "usage: %prog [options] < fastaInput"
optparser = OptionParser(usage)
optparser.add_option("-b", "--binary", dest="bin",
                     help="Path to genscan binary",
                     default="genscan")
optparser.add_option("-s", "--seglength", type="int", dest="segLength",
                     help="break large sequences into SEGLENGTH size pieces",
                     default=2500000)
optparser.add_option("-p", "--paramfile", dest="paramFilename",
                     help="use PARAMFILE for genscan parameters",
                     default="/usr/local/apps/genscan/HumanIso.smat",
                     metavar="PARAMFILE")
(options, args) = optparser.parse_args()
if len(args) != 0:
    optparser.error("incorrect number of arguments")

def isAllNs(seq):
    return seq.upper() == ("N" * len(seq))

def runGenscan(rec, bin, segLength, paramFilename):
    seqFilename = os.tmpnam()

    cmd = "%(bin)s %(paramFilename)s %(seqFilename)s | genscan2gtf -s%(segLength)d" % vars()

    gffRecs = []
    subrec = FASTA.Record()
    for i in range(((len(rec.sequence) - 1)/ segLength) + 1):
        subrec.title = "%s.%d" % (rec.title, i)
        subrec.sequence = rec.sequence[i * segLength: (i + 1) * segLength]
        if isAllNs(subrec.sequence):
            continue
        seqFile = file(seqFilename, 'w')
        seqFile.write(str(subrec))
        seqFile.close()
        gffRecs.extend([gffRec for gffRec in GFF.Iterator(os.popen(cmd))])

    try:
        os.remove(seqFilename)
    except OSError, e:
        sys.stderr.write("Could not delete temporary file %s: %s" % \
                         (seqFilename, str(e)))
    return gffRecs

for rec in FASTA.Iterator(sys.stdin):
    gffRecs = runGenscan(rec, options.bin, options.segLength, options.paramFilename)
    for gffRec in gffRecs:
        print gffRec
