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
import GFF
from optparse import OptionParser

usage = "usage: %prog bad.gff < gffInput"
optparser = OptionParser(usage)
optparser.add_option("-a", "--attr", dest="attr", default="name")
(options, args) = optparser.parse_args()
if len(args) != 1:
    optparser.error("Invalid number of arguments")
badGFFFilename = args[0]

bad = {}
attr = options.attr
for rec in GFF.Iterator(file(badGFFFilename)):
    if attr in rec.attributes:
        bad[rec.attributes[attr]] = None

for rec in GFF.Iterator(sys.stdin):
    if attr not in rec.attributes or rec.attributes[attr] not in bad:
        print rec

        
