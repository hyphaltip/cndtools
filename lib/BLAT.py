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

class Record:
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

    def __str__(self):
        return '\t'.join(map(str, (self.queryID,
                                   self.subjectID,
                                   self.pctIdentity,
                                   self.alignLength,
                                   self.mismatches,
                                   self.gaps,
                                   self.queryStart,
                                   self.queryEnd,
                                   self.subjectStart,
                                   self.eValue,
                                   self.score)))

class Iterator:
    def __init__(self, handle):
        self.handle = handle

    def __iter__(self):
        return self

    def next(self):
        while 1:
            line = self.handle.readline()
            if not line:
                raise StopIteration
            fields = line[:-1].split('\t')

            r = Record()
            r.queryID, r.subjectID = fields[0:2]
        
            r.pctIdentity, r.eValue = map(float, (fields[2], fields[10]))
            (r.alignLength, r.mismatches, r.gaps, r.queryStart, r.queryEnd,
             r.subjectStart, r.subjectEnd) = map(int, fields[3:10])
            r.score = int(float(fields[11]))
        
            return r
