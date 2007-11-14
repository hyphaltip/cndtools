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
from optparse import OptionParser

import FASTA
import ZFF

usage = "usage: %prog [options] hmmFile < fastaInput"
optparser = OptionParser(usage)
optparser.add_option("-b", "--binary", dest="bin",
                     help="Path to SNAP binary",
                     default="snap")
optparser.add_option("-s", "--seglength", type="int", dest="segLength",
                     help="break large sequences into SEGLENGTH size pieces",
                     default=10000000)
optparser.add_option("--lcmask", action="store_true", dest="lcmask", default=0)
(options, args) = optparser.parse_args()
if len(args) != 1:
    optparser.error("incorrect number of arguments")
hmmFile = args[0]

def runSnap(rec, bin, segLength, paramFilename):
    snapOpts = ""
    if options.lcmask:
        snapOpts += " -lcmask"
    cmd = "%(bin)s %(snapOpts)s %(paramFilename)s /dev/stdin"

    subrec = FASTA.Record()
    for i in range(((len(rec.sequence) - 1)/ segLength) + 1):
        subrec.title = "%s.%d" % (rec.title, i)
        subrec.sequence = rec.sequence[i * segLength: (i + 1) * segLength]
        print >>sys.stderr, subrec.title
        snapIn, snapOut = os.popen2(cmd % vars())
        snapIn.write(str(subrec))
        snapIn.close()
        for zffRec in ZFF.Iterator(snapOut):
            splitIndex = zffRec.seqname.rfind(".")
            segNum = int(zffRec.seqname[splitIndex + 1:])
            offset = segNum * segLength
            zffRec.seqname = zffRec.seqname[:splitIndex]
            for feature in zffRec.features:
                feature.start += offset
                feature.end += offset
            print zffRec

for rec in FASTA.Iterator(sys.stdin):
    runSnap(rec, options.bin, options.segLength, hmmFile)
