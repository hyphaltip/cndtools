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
import os
import gc
from optparse import OptionParser

usage = "usage: %prog [options] > constraints"
optparser = OptionParser(usage)
optparser.add_option("--noncoding", dest="outputNonCoding",
                     action="store_true",
                     help="Output constraints for non coding hits",
                     default=0, metavar="PATH")
optparser.add_option("-i", "--input-dir", dest="inputDir",
                     help="Directory containing hit and BLAT files",
                     default=".", metavar="PATH")
optparser.add_option("-m", "--map-dir", dest="mapDir",
                     help="Directory containing map related files",
                     default=".", metavar="PATH")
(options, args) = optparser.parse_args()

if len(args) != 0:
    optparser.error("Invalid number of arguments")

class BlatHit:
    # Fields:
    # queryID
    # subjectID
    # pctIdentity
    # alignLength
    # mismatches
    # gaps
    # queryStart
    # queryEnd
    # subjectStart
    # subjectEnd
    # eValue
    # score
    pass

class BlatHitFileParser:
    def __init__(self, handle, phits):
        self.handle = handle
        self.phits = phits

    def __iter__(self):
        return self

    def next(self):
        while 1:
            line = self.handle.readline()
            if not line:
                raise StopIteration
            fields = line[:-1].split('\t')

            h = BlatHit()
            h.queryID, h.subjectID = map(int, fields[0:2])

            try:
                if self.phits[h.queryID].anchor2 != h.subjectID:
                    continue
            except:
                continue
        
            h.pctIdentity, h.eValue = map(float, (fields[2], fields[10]))
            (h.alignLength, h.mismatches, h.gaps, h.queryStart, h.queryEnd,
             h.subjectStart, h.subjectEnd) = map(long, fields[3:10])
            h.score = fields[11]
        
            return h

class PairwiseHit:
    def __init__(self, run_id, genome1, anchor1, genome2, anchor2):
        if genome1 > genome2:
            genome1, genome2 = genome2, genome1
            anchor1, anchor2 = anchor2, anchor1
        self.run_id = int(run_id)
        self.genome1 = genome1
        self.genome2 = genome2
        self.anchor1 = anchor1
        self.anchor2 = anchor2

class PairwiseHitFileParser:
    def __init__(self, handle):
        self.handle = handle

    def __iter__(self):
        return self

    def next(self):
        line = self.handle.readline()
        if not line:
            raise StopIteration
        return PairwiseHit(*line.split())

class ProteinAnchor:
    # Fields:
    # id
    # chrom
    # strand
    # start
    # end
    # isCoding
    pass

class ProteinAnchorFileParser:
    def __init__(self, handle):
        self.handle = handle

    def __iter__(self):
        return self

    def next(self):
        line = self.handle.readline()
        if not line:
            raise StopIteration
        fields = line[:-1].split('\t')
        p = ProteinAnchor()
        p.id = int(fields[0])
        p.chrom, p.strand = fields[1:3]
        p.start, p.end = map(int, fields[3:5])
        p.isCoding = int(fields[5])
        return p

def readPairwiseHitFile(handle):
    d = {}
    for ph in PairwiseHitFileParser(handle):
        hits = d.setdefault(ph.genome1, {}).setdefault(ph.genome2, {})
        assert(ph.anchor1 not in hits)
        hits[ph.anchor1] = ph
    return d

def readProteinAnchorFile(handle):
    d = {}
    for anchor in ProteinAnchorFileParser(handle):
        d[anchor.id] = anchor
    return d

def hitsConsistent(hit1, hit2):
    return (hit1[0] < hit2[0]) == (hit1[1] < hit2[1])

def isConsistentWith(hit, hits):
    for hit2 in hits:
        if not hitsConsistent(hit, hit2):
            return 0
    return 1

def intervalsOverlap((start1, end1), (start2, end2)):
    return start1 <= end2 and end1 >= start2

def overlapsWithAny(i, intervals):
    for i2 in intervals:
        if intervalsOverlap(i, i2):
            return 1
    return 0

def getHits(handle, hits, anchors):
    queryIntervalsMap = {}
    subjectIntervalsMap = {}
    for bhit in BlatHitFileParser(handle, hits):
        phit = hits[bhit.queryID]

        query = anchors[phit.genome1][bhit.queryID]
        subject = anchors[phit.genome2][bhit.subjectID]

        queryIntervals = queryIntervalsMap.setdefault(bhit.queryID, [])
        subjectIntervals = subjectIntervalsMap.setdefault(bhit.subjectID, [])
            
        # Construct intervals used in this alignment
        queryInterval = (bhit.queryStart, bhit.queryEnd)
        subjectInterval = (bhit.subjectStart, bhit.subjectEnd)

        # Check that this alignment does not overlap with a
        # previous alignment
        if (overlapsWithAny(queryInterval, queryIntervals) or
            overlapsWithAny(subjectInterval, subjectIntervals)):
            continue

        # Check that this alignment is consistent with the previous
        # alignments
        if (not isConsistentWith((queryInterval, subjectInterval),
                                 zip(queryIntervals, subjectIntervals))):
            continue

        queryIntervals.append(queryInterval)
        subjectIntervals.append(subjectInterval)

        if query.strand == '+':
            queryStart = (bhit.queryStart - 1) * 3 + query.start
            queryEnd = bhit.queryEnd * 3 + query.start
        else:
            queryStart = query.end - bhit.queryEnd * 3
            queryEnd = query.end - (bhit.queryStart - 1) * 3
            
        if subject.strand == '+':
            subjectStart = (bhit.subjectStart - 1) * 3 + subject.start
            subjectEnd = bhit.subjectEnd * 3 + subject.start
        else:
            subjectStart = subject.end - bhit.subjectEnd * 3
            subjectEnd = subject.end - (bhit.subjectStart - 1) * 3
            
        print '\t'.join(map(str, (phit.run_id,
                                  1,
                                  phit.genome1,
                                  query.chrom,
                                  queryStart,
                                  queryEnd,
                                  query.strand,
                                  phit.genome2,
                                  subject.chrom,
                                  subjectStart,
                                  subjectEnd,
                                  subject.strand)))

def getNonCodingHits(hits, anchors):
    for phit in hits.values():
        a1 = anchors[phit.genome1][phit.anchor1]
        a2 = anchors[phit.genome2][phit.anchor2]
        if a1.isCoding or a2.isCoding: continue
        print '\t'.join(map(str, (phit.run_id,
                                  0,
                                  phit.genome1,
                                  a1.chrom, a1.start, a1.end, a1.strand,
                                  phit.genome2,
                                  a2.chrom, a2.start, a2.end, a2.strand)))

genomes = file(os.path.join(options.mapDir, "genomes")).read().split()

sys.stderr.write("Reading pairwise hit file...")
phitFilename = os.path.join(options.mapDir, "pairwisehits")
phits = readPairwiseHitFile(open(phitFilename))
sys.stderr.write("done.\n")
gc.collect()


sys.stderr.write("Reading anchor files...")
anchors = {}
for genome in genomes:
    sys.stderr.write(" " + genome)

    # Test for existance of draft genome anchor file
    draftAnchorFilename = os.path.join(options.mapDir, genome + ".anchors")
    if os.path.exists(draftAnchorFilename):
        anchorFilename = draftAnchorFilename
    else:
        anchorFilename = os.path.join(options.inputDir, genome + ".anchors")
    anchors[genome] = readProteinAnchorFile(file(anchorFilename))
    gc.collect()
sys.stderr.write(" done.\n")

sys.stderr.write("Reading blat files...\n")
for (genome1, genomeHits) in phits.items():
    for (genome2, hits) in genomeHits.items():
        blatFileName = "%s-%s.blat" % (genome1, genome2)
        blatFile = file(os.path.join(options.inputDir, blatFileName))
        sys.stderr.write(blatFileName + "\n")
        getHits(blatFile, hits, anchors)
        if options.outputNonCoding:
            getNonCodingHits(hits, anchors)
