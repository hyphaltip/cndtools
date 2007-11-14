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


from Bio import GenBank
from Bio.WWW import RequestLimiter
import sys
import math

if len(sys.argv) != 5:
    print "Usage: %s accessionPrefix startNum endNum digits" % sys.argv[0]
    sys.exit(1)

# Extract command line arguments
prefix = sys.argv[1]
start, end, digits = map(long, sys.argv[2:])

# Receiver functor - prints out records passed to it to HANDLE
class RecordReceiver:
    def __init__(self, handle):
        self.handle = handle
    def __call__(self, id, rec):
        self.handle.write(rec)

# Functor that deals with bad records - prints an error message to HANDLE
class BadRecordReceiver:
    def __init__(self, handle):
        self.handle = handle
        self.badIDs = []
    def __call__(self, badID):
        self.badIDs.append(badID)
        self.handle.write("Bad ID: %s\n" % badID)

# Form pattern for accession strings
accessionPat = prefix + "%%0%dd" % digits

batchSize = 500

for curr in xrange(start, end + 1, batchSize):
    # Generate accession strings for this batch
    ids = [accessionPat % num for num in range(curr,
                                               min(curr + batchSize, end + 1))]
    GenBank.download_many(ids,
                          RecordReceiver(sys.stdout),
                          database='nucleotide',
                          broken_fn=BadRecordReceiver(sys.stderr),
                          faildelay=5.0)

