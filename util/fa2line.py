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


import FASTA
import sys

recs = FASTA.readFasta(sys.stdin)

maxLen = max([len(rec.sequence) for rec in recs])

titleWidth = 10
seqWidth = 60

start = 0
while start < maxLen:
    for rec in recs:
        if start == 0:
            print "%s%s" % (rec.title.ljust(titleWidth)[0:titleWidth],
                            rec.sequence[start : start + seqWidth])
        else:
            print "%s%s" % (' ' * titleWidth,
                            rec.sequence[start : start + seqWidth])
    start += seqWidth
    print
