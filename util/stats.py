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
import stats
from optparse import OptionParser

usage = "usage: %prog [options] < numericalInput"
parser = OptionParser(usage)
parser.add_option("-f", "--float", dest="useFloatValues", action="store_true",
                  help="input consists of floating point numbers", default=0)
parser.add_option("--mode", dest="showMode", action="store_true",
                  help="show mode of input values", default=0)
(options, args) = parser.parse_args()
if len(args) != 0:
    parser.error("incorrect number of arguments")

try:
    if options.useFloatValues:
        d = map(float, sys.stdin.xreadlines())
    else:
        d = map(long, sys.stdin.xreadlines())
except ValueError, err:
    sys.stderr.write("Bad datum: %s\n" % str(err))
    sys.exit(1)

if len(d) == 0:
    sys.stderr.write("No data given\n")
    sys.exit(1)
        
d.sort()

print "           N =", len(d)
print "         SUM =", stats.sum(d)
print "         MIN =", min(d)
print "1ST-QUARTILE =", stats.firstquartile(d)
print "      MEDIAN =", stats.median(d)
print "3RD-QUARTILE =", stats.thirdquartile(d)
print "         MAX =", max(d)
print "        MEAN =", stats.mean(d)

if d[0] < 0:
    print "         N50 = NA"
else:
    print "         N50 =", stats.n50(d)

if options.showMode:
    modeinfo = stats.mode(d)
    print "     MODE(S) =", ','.join(map(str, modeinfo[0])), "(%d)" % modeinfo[1]
