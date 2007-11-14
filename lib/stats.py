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

import operator

def sum(L): return reduce(operator.add, L, 0)

def product(L): return reduce(operator.mul, L, 1)

def mean(L): return float(sum(L)) / len(L)

def midval(L, low, high):
    numelts = high - low
    if numelts % 2 == 1:
        return L[low + numelts / 2]
    else:
        return mean(L[low - 1 + numelts / 2 : low + 1 + numelts / 2])

def median(L):
    return midval(L, 0, len(L))

def n50(L):
    fiftypct = sum(L) / 2
    partialsum = 0
    i = len(L)
    while partialsum < fiftypct:
        i -= 1
        partialsum += L[i]
    return L[i]

def firstquartile(L):
    return midval(L, 0, len(L) / 2)

def thirdquartile(L):
    return midval(L, len(L) / 2, len(L))

def counts(L):
    cts = {}
    for x in L:
        cts[x] = cts.setdefault(x, 0) + 1
    return cts
    
def mode(L):
    cts = counts(L)
    countpairs = zip(cts.values(), cts.keys())
    countpairs.sort()
    maxcount = countpairs[-1][0]
    return ([key for ct, key in countpairs if ct == maxcount], maxcount)
