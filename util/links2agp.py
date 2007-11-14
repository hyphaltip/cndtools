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


from optparse import OptionParser
import sys

import AGP

usage = "usage: %prog < linksInput > agpOutput"
optparser = OptionParser(usage)
optparser.add_option("-g", "--gap-length", type="long",
                     dest="gapLength", metavar="LENGTH", default=25,
                     help="Gap length to use when estimated gap length is negative")
optparser.add_option("--scaffold-prefix", dest="scaffoldPrefix",
                     metavar="STRING", default="scaffold_",
                     help="Prefix for scaffold names")
optparser.add_option("--contig-prefix", dest="contigPrefix",
                     metavar="STRING", default="contig_",
                     help="Prefix for contig names")

(options, args) = optparser.parse_args()
if len(args) != 0:
    optparser.error("Invalid number of arguments")

scaffoldContigs = {}

for line in sys.stdin:
    if line.startswith('#'): continue

    fields = line.split()
    numContigInScaffold = int(fields[3])
    contigID = fields[4]
    contigLength = long(fields[5])
    gapBeforeContig = long(fields[6])
    scaffoldContigs.setdefault(fields[0], []).append((numContigInScaffold,
                                                      contigID,
                                                      contigLength,
                                                      gapBeforeContig))


rec = AGP.Record()
rec.recNum = 1

# Fields that are the same for all gaps
rec.gapKind = "contig"
rec.isBridged = 1

# Fields that are the same for all contig records
rec.sourceStart = 1
rec.sourceOrientation = '+'

scaffoldContigPairs = scaffoldContigs.items()
scaffoldContigPairs.sort()

for scaffoldID, contigs in scaffoldContigPairs:
    contigs.sort()
    rec.chrom = options.scaffoldPrefix + scaffoldID
    rec.contigStart = 1
    for num, contigID, contigLength, gapBeforeContig in contigs:
        if gapBeforeContig < 0:
            gapBeforeContig = options.gapLength
        if gapBeforeContig > 0:
            rec.contigEnd = rec.contigStart + gapBeforeContig - 1
            rec.type = 'N'
            rec.gapLength = gapBeforeContig
            print str(rec)
            
            rec.recNum += 1
            rec.contigStart = rec.contigEnd + 1

        rec.contigEnd = rec.contigStart + contigLength - 1
        rec.type = 'D'
        rec.sourceAccession = options.contigPrefix + contigID
        rec.sourceEnd = contigLength
        print str(rec)

        rec.recNum += 1
        rec.contigStart = rec.contigEnd + 1


                                                      
                                                     
