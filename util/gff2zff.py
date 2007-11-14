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


from optparse import OptionParser
import sys

import GFF
import ZFF

usage = "usage: %prog < input > output"
optparser = OptionParser(usage)
optparser.add_option("-a", "--attr", dest="attr", default="gene_id")
(options, args) = optparser.parse_args()
if len(args) != 0:
    optparser.error("Invalid number of arguments")

genes = {}
geneFeatures = ('start_codon', 'stop_codon', 'CDS')
nameAttr = options.attr

for rec in GFF.Iterator(sys.stdin):
    if rec.feature in geneFeatures and nameAttr in rec.attributes:
        chromGenes = genes.setdefault(rec.seqname, {})
        gene = chromGenes.setdefault(rec.attributes[nameAttr], [])
        gene.append(rec)

pairs = genes.items()
pairs.sort()

for seqname, seqGenes in pairs:
    print ">%s" % seqname
    for geneName, recs in seqGenes.items():
        features = [rec.feature for rec in recs]
        if 'start_codon' not in features or 'stop_codon' not in features:
            continue
        cdsRecs = [rec for rec in recs if rec.feature == "CDS"]
        assert(len(cdsRecs) >= 1)
        strand = cdsRecs[0].strand

        if strand == '+':
            cds = [(rec.start, rec.end) for rec in cdsRecs]
        else:
            cds = [(rec.end, rec.start) for rec in cdsRecs]

        cds.sort()

        if len(cds) == 1:
            print ZFF.Record("Esngl", cds[0][0], cds[0][1], geneName)
        else:
            if strand == '+':
                print ZFF.Record("Einit", cds[0][0], cds[0][1], geneName)
            else:
                print ZFF.Record("Eterm", cds[0][0], cds[0][1], geneName)
                
            for begin, end in cds[1:-1]:
                print ZFF.Record("Exon", begin, end, geneName)

            if strand == '+':
                print ZFF.Record("Eterm", cds[-1][0], cds[-1][1], geneName)
            else:
                print ZFF.Record("Einit", cds[-1][0], cds[-1][1], geneName)                
