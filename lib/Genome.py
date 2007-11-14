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

import os
import string
import FASTA
import DB_MySQL

genomeDir = os.environ["GENOME_PATH"]

def setGenomeDir(d):
    global genomeDir
    genomeDir = d

def getDBSeq(db, chrom, start, end,
             strand='+', masked=0, hardMasked=0):
    options = []
    if hardMasked:
        options.append('-h')
    elif not masked:
        options.append('-u')
    options = ' '.join(options)

    cmd = ("sdbExport %(options)s %(db)s %(chrom)s " + 
           "%(start)d %(end)d %(strand)s") % vars()
    file = os.popen(cmd, 'r')
    return FASTA.Iterator(file).next().sequence

def getGenomicSeq(genome, chrom, start, end,
                  strand='+', masked=0, hardMasked=0):
    dbFilename = os.path.join(genomeDir, genome + ".sdb")
    return getDBSeq(dbFilename, chrom, start, end, strand, masked, hardMasked)

def getChromLen(genome, chrom):
    nibFilename = _getChromFilename(genome, chrom)
    cmd = ("nibLen %(nibFilename)s" % vars())
    return int(os.popen(cmd, 'r').readline())

def importSeq(genome, id, seq):
    DB = DB_MySQL.MySQLDB(db=genome, passwd="stupid")
    DB.insertRec("seq", (id, seq))

