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


import BLAT
import sys
import stats
from optparse import OptionParser

usage = "usage: %prog [options] < blatInput > hitsOutput"
optparser = OptionParser(usage)
optparser.add_option("-i", "--identity", action="store_true", dest="identity",
                     help="remove hits with identical subject and query",
                     default=0)
optparser.add_option("--no-combine", action="store_false", dest="combine",
                     help="do not combine multiple hits between the same subject and query",
                     default=1)
(options, args) = optparser.parse_args()
if len(args) != 0:
    optparser.error("incorrect number of arguments")


def intervalsOverlap((start1, end1), (start2, end2)):
    return start1 <= end2 and end1 >= start2

def overlapsWithAny(i, intervals):
    for i2 in intervals:
        if intervalsOverlap(i, i2):
            return 1
    return 0

def tabjoin(seq):
    return '\t'.join(map(str, seq))

if not options.combine:
    if not options.identity:
        for hit in BLAT.Iterator(sys.stdin):
            print tabjoin((hit.queryID, hit.subjectID, hit.score, hit.eValue))
    else:
        for hit in BLAT.Iterator(sys.stdin):
            if hit.queryID != hit.subjectID:
                print tabjoin((hit.queryID, hit.subjectID, hit.score, hit.eValue))                
else:
    lastQuery = None

    for hit in BLAT.Iterator(sys.stdin):
        q = hit.queryID
        s = hit.subjectID
        qInterval = (hit.queryStart, hit.queryEnd)
        sInterval = (hit.subjectStart, hit.subjectEnd)

        if options.identity and q == s:
            continue
        elif q != lastQuery:
            if lastQuery is not None:
                for subject in scores:
                    score = stats.sum(scores[subject])
                    eValue = stats.product(eValues[subject])
		    if eValue < 1e-50:
			eValue = 0
                    print tabjoin((lastQuery, subject, score, "%.0e" % eValue))

            lastQuery = q
            queryIntervals = {s: [qInterval]}
            subjectIntervals = {s: [sInterval]}
            scores = {s: [hit.score]}
            eValues = {s: [hit.eValue]}

        elif ((not overlapsWithAny(qInterval,
                                   queryIntervals.setdefault(s, []))) and
              (not overlapsWithAny(sInterval,
                                   subjectIntervals.setdefault(s, [])))):
            queryIntervals[s].append(qInterval)
            subjectIntervals[s].append(sInterval)
            scores.setdefault(s, []).append(hit.score)
            eValues.setdefault(s, []).append(hit.eValue)
    
    if lastQuery is not None:
        for subject in scores:
            score = stats.sum(scores[subject])
            eValue = stats.product(eValues[subject])
	    if eValue < 1e-50:
		eValue = 0
            print tabjoin((lastQuery, subject, score, "%.0e" % eValue))
