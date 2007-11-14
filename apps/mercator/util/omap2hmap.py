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
import TabDelimited
from optparse import OptionParser

usage = "usage: %prog [options] genomeFile < orthologyMap > homologyMap"
optparser = OptionParser(usage)
(options, args) = optparser.parse_args()
if len(args) != 1:
    optparser.error("incorrect number of arguments")

genomeFilename = args[0]

genomes = file(genomeFilename).read().split()

for fields in TabDelimited.Iterator(sys.stdin):
    newfields = [fields[0]]
    for i in xrange(len(genomes)):
        chrom, start, end, strand = fields[i * 4 + 1: i * 4 + 5]
        if chrom != "NA":
            newfields.extend((genomes[i], chrom, start, end, strand))
    print '\t'.join(newfields)
