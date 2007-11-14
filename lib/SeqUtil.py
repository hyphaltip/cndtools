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
import math
import types

def isTupleOrList(x):
    return (type(x) in (types.TupleType, types.ListType))

def flatten(seq):
    flattenedSeq = []
    for elt in seq:
        if isTupleOrList(elt):
            flattenedSeq.extend(flatten(elt))
        else:
            flattenedSeq.append(elt)
    return flattenedSeq

def sum(seq):
    """Returns the sum of the elements in a sequence"""
    return reduce(operator.add, seq, 0)

def product(seq):
    """Returns the product of the elements in a sequence"""
    return reduce(operator.mul, seq, 1)

def mean(seq):
    """Returns the arithmetic mean of the elements in a sequence"""
    return float(sum(seq)) / len(seq)

def var(seq):
    """Returns the variance of the elements in a sequence"""
    mu = mean(seq)
    sigmaSquared = 0
    for elt in seq:
        sigmaSquared += (elt - mu) * (elt - mu)
    return float(sigmaSquared) / len(seq)

def std(seq):
    """Returns the standard deviation of the elements in a sequence"""
    return math.sqrt(var(seq))

def median(seq):
    """Returns the median of the elements in a sequence"""
    sortedSeq = seq[:]
    sortedSeq.sort()
    middle, odd = divmod(len(seq), 2)
    if odd:
        return sortedSeq[len(seq) / 2]
    else:
        return mean(sortedSeq[middle - 1: middle + 1])

def mode(seq):
    """Returns the mode of the elements in a sequence"""
    uniqueElts = {}
    for elt in seq:
        uniqueElts[elt] = uniqueElts.get(elt, 0) + 1
    maxCount = max(uniqueElts.values())
    return [elt for (elt, count) in uniqueElts.items() if count == maxCount]

def longestIncrSubseq(seq):
    """Returns a longest increasing subsequence of a sequence"""

    # Check for empty sequence
    if len(seq) == 0:
        return []

    # Create greedy cover of decreasing subsequences
    cover = []
    for elt in seq:
        # Find leftmost decreasing subseq to append elt to
        for subSeq in cover:
            if (elt <= subSeq[-1]):
                subSeq.append(elt)
                break
        else:
            cover.append([elt])
    
    # Create increasing subsequence
    cover.reverse()
    x = cover[0][0]
    lis = [x]
    for subSeq in cover[1:]:
        # Pick an element from each decreasing subsequence
        for elt in subSeq:
            # Find first number that is smaller than x
            if elt < x:
                lis.append(elt)
                x = elt
                break
    lis.reverse()
    return lis

