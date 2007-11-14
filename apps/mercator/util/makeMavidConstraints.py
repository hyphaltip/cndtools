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

usage = "usage: %prog [options] multi.fa workdir"
optparser = OptionParser(usage)
optparser.add_option("--genscan-params", dest="genscanParamFile",
                     help="Full path to parameter file to be used by genscan",
                     default="HumanIso.smat", metavar="FILEPATH")
(options, args) = optparser.parse_args()
if len(args) != 2:
    optparser.error("incorrect number of arguments")

multiFastaFilename, workdir = args

def firstWord(s):
    words = s.strip().split()
    return words[0]

def readAnchors(strm):
    anchors = {}
    for line in strm:
        fields = line.split()
        anchors[fields[0]] = fields[3:5]
    return anchors

# Test for existence of genscan parameter file
if not os.path.exists(options.genscanParamFile):
    sys.stderr.write("Error: Genscan parameter file %s does not exist\n" %
                     options.genscanParamFile)
    sys.exit(1)

# Read in the sequences from the multi-fasta file
sys.stderr.write("Reading in multi-fasta file...")
multiFastaFile = file(multiFastaFilename)
fastaRecs = FASTA.readFasta(multiFastaFile)
multiFastaFile.close()
sys.stderr.write("done\n")

# Make protein anchors for each sequence
for rec in fastaRecs:
    rec.title = firstWord(rec.title)
    chromFile = file(os.path.join(workdir, rec.title + ".chroms"), 'w')
    chromFile.write("%s\t%d\n" % (rec.title, len(rec.sequence)))
    chromFile.close()

    sys.stderr.write("Writing single-fasta file...")
    fastaFilename = os.path.join(workdir, rec.title + ".fa")
    fastaFile = file(fastaFilename, 'w')
    fastaFile.write(str(rec))
    fastaFile.close()
    sys.stderr.write("done\n")

    sys.stderr.write("Running genscan...")
    gffFilename = os.path.join(workdir, rec.title + ".gff")
    os.system("runGenscan -p%s < %s > %s" %
              (options.genscanParamFile, fastaFilename, gffFilename))
    sys.stderr.write("done\n")

    sys.stderr.write("Making protein anchors...")
    anchorFilename = os.path.join(workdir, rec.title + ".anchors")
    os.system("gff2anchors < %s > %s" % (gffFilename, anchorFilename))
    sys.stderr.write("done\n")

    sys.stderr.write("Extracting protein sequences for anchors...")
    proteinFilename = os.path.join(workdir, rec.title + ".proteins.fa")
    os.system("anchors2fa --fasta %s < %s > %s" %
              (fastaFilename, anchorFilename, proteinFilename))
    sys.stderr.write("done\n")
    
# BLAT protein anchors pairwise
sys.stderr.write("BLATing anchors pairwise...\n")
for rec1 in fastaRecs:
    for rec2 in fastaRecs:
        if rec1.title > rec2.title:
            continue
        pairName = rec1.title + "-" + rec2.title
        sys.stderr.write(pairName + "\n")
        protein1 = os.path.join(workdir, "%s.proteins.fa" % rec1.title)
        protein2 = os.path.join(workdir, "%s.proteins.fa" % rec2.title)
        blatOutput = os.path.join(workdir, "%s.blat" % pairName)
        os.system("blat -prot -out=blast8 %s %s %s" %
                  (protein2, protein1, blatOutput))
        hitOutput = os.path.join(workdir, "%s.hits" % pairName)
        os.system("blat2hits < %s > %s" % (blatOutput, hitOutput))

# Make map
sys.stderr.write("Making map...\n")
cmd = "mercator -i %s -o %s %s" % \
      (workdir, workdir, ' '.join([rec.title for rec in fastaRecs]))
sys.stderr.write(cmd + '\n')
os.system(cmd)

# Make constraints file
sys.stderr.write("Making constraints file...\n")
constraintsFilename = os.path.join(workdir, "constraints")
cmd = "phits2constraints %s %s > %s" % \
      (workdir, workdir, constraintsFilename)
os.system(cmd)

# Make MAVID constraints
sys.stderr.write("Making MAVID constraints...\n")
mavidConstraintsFile = file(os.path.join(workdir, "mavidconstraints"), 'w')
for line in file(constraintsFilename):
    fields = line.split()
    if fields[3] == fields[8]:
        print >>mavidConstraintsFile, \
              '\t'.join([fields[1], fields[4], fields[5],
                         fields[6], fields[9], fields[10]])

# Make MAVID runs
sys.stderr.write("Making MAVID runs...\n")
runsFile = file(os.path.join(workdir, "runs"))
mavidRunsFile = file(os.path.join(workdir, "mavidruns"), 'w')
genomes = [rec.title for rec in fastaRecs]
anchors = [readAnchors(file(os.path.join(workdir, genomes[i] + ".anchors")))
           for i in xrange(len(genomes))]
for line in runsFile:
    fields = line.split()
    coords = []
    for genome, anchor in zip(xrange(len(genomes)), fields):
        if anchor == "NA":
            coords.extend(["-1", "-1"])
        else:
            coords.extend(anchors[genome][anchor])
    print >>mavidRunsFile, '\t'.join(coords)
