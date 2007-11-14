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
from optparse import OptionParser

import FASTA
import GFF

usage = "usage: %prog [options] < fastaInput > gffOutput"
optparser = OptionParser(usage)
optparser.add_option("-s", "--seglength", type="int", dest="segLength",
                     help="break large sequences into SEGLENGTH size pieces",
                     default=100000000)
optparser.add_option("-p", "--paramfile", dest="paramFilename",
                     help="use PARAMFILE for parameters",
                     default="/usr/local/apps/augustus/config/human_singlestrand_partial.cfg",
                     metavar="PARAMFILE")
(options, args) = optparser.parse_args()
if len(args) != 0:
    optparser.error("incorrect number of arguments")

def runAugustus(rec,
                param="/usr/local/apps/augustus/config/human_singlestrand_partial.cfg",
                segLength=500000,
                bin="/usr/local/apps/augustus/bin/augustus",
                options=None):
    if not options:
        options = []

    seqFilename = os.tmpnam()

    options.append("--strand=both")
    options.append("--protein=off")
    options.append("--start=off")
    options.append("--stop=off")
    options.append("--queryfile=" + seqFilename)
    options.append("--/NAMGene/maxDNAPieceSize=" + str(segLength))
    
    optString = ' '.join(options)

    cmd = "%(bin)s %(param)s %(optString)s" % vars()

    gffRecs = []

    subrec = FASTA.Record()
    subrec.title = rec.title
    for i in range(((len(rec.sequence) - 1)/ segLength) + 1):
        subrec.sequence = rec.sequence[i * segLength: (i + 1) * segLength]
        seqFile = file(seqFilename, 'w')
        seqFile.write(str(subrec))
        seqFile.close()

        for line in os.popen(cmd):
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')

            cdsRec = GFF.Record(seqname=rec.title,
                                source="geneid",
                                feature="CDS",
                                start=int(fields[3]) + i * segLength,
                                end=int(fields[4]) + i * segLength,
                                score=float(fields[5]),
                                strand=fields[6],
                                frame=fields[7],
                                attributes={"gene_id": [fields[8]],
                                            "transcript_id": [fields[8] + ".1"]})

            exonType = fields[2]

            if exonType in ["First", "Single"]:
                startCodonRec = cdsRec.copy()
                startCodonRec.feature = "start_codon"
                startCodonRec.score = None
                startCodonRec.frame = None
                if cdsRec.strand == '+':
                    startCodonRec.end = startCodonRec.start + 2
                else:
                    startCodonRec.start = startCodonRec.end - 2
                gffRecs.append(startCodonRec)

            exonRec = cdsRec.copy()
            exonRec.feature = "exon"
            exonRec.frame = None
            gffRecs.append(exonRec)

            gffRecs.append(cdsRec)

            if exonType in ["Terminal", "Single"]:
                stopCodonRec = cdsRec.copy()
                stopCodonRec.feature = "stop_codon"
                stopCodonRec.score = None
                stopCodonRec.frame = None
                if cdsRec.strand == '+':
                    stopCodonRec.start = stopCodonRec.end - 2
                else:
                    stopCodonRec.end = stopCodonRec.start + 2
                gffRecs.append(stopCodonRec)
                
    try:
        os.remove(seqFilename)
    except OSError, e:
        sys.stderr.write("Could not delete temporary file %s: %s" % \
                         (seqFilename, str(e)))

    return gffRecs         

for rec in FASTA.Iterator(sys.stdin):
    print >>sys.stderr, rec.title
    gffRecs = runGeneid(rec,
                        param=options.paramFilename,
                        segLength=options.segLength)
    for gffRec in gffRecs:
        print gffRec
