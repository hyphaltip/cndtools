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
import sys
import FASTA

dnaAmbiguities = {
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "N": "GATC",
    }

def unambiguousChar(c, ambiguousChars=dnaAmbiguities):
    if c in ambiguousChars:
        return random.choice(ambiguousChars[c])
    else:
        return c

def unambiguousSeq(seq):
    return ''.join(map(unambiguousChar, seq))

for rec in FASTA.Iterator(sys.stdin):
    rec.sequence = unambiguousSeq(rec.sequence)
    print rec

