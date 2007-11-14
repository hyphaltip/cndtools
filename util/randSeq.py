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


import random
import FASTA
import sys

if len(sys.argv) != 3:
    print >> sys.stderr, "usage: %s title length" % sys.argv[0]
    sys.exit(1)

title = sys.argv[1]
length = int(sys.argv[2])

DNA = "ACGT"
seq = [None] * length
for i in xrange(length):
    seq[i] = random.choice(DNA)

rec = FASTA.Record()
rec.title = title
rec.sequence = ''.join(seq)
print rec

