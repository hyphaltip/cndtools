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
    """Holds information from a FASTA record.

    Members:
    title       Title line ('>' character not included).
    sequence    The sequence.
    
    """
    def __init__(self, title="", sequence="", colwidth=60):
        """__init__(self, colwidth=60)

        Create a new Record.  colwidth specifies the number of residues
        to put on each line when generating FASTA format.

        """
        self.title = title
        self.sequence = sequence
        self._colwidth = colwidth
        
    def __str__(self):
        s = []
        s.append('>%s' % self.title)
        i = 0
        while i < len(self.sequence):
            s.append(self.sequence[i:i+self._colwidth])
            i = i + self._colwidth
        return '\n'.join(s)

def readFasta(strm):
    txt = strm.read()
    entries = []
    for entry in txt.split('>')[1:]:
        rec = Record()
        rec.title, rec.sequence = entry.split('\n', 1)
        rec.sequence = rec.sequence.replace('\n', '').replace(' ', '').upper()
        entries.append(rec)
    return entries

class Iterator:
    def __init__(self, stream, bufsize=(100 * pow(2, 20))):
        self.stream = stream
        self.bufsize = bufsize
        self.pos = -1
        while self.pos == -1:
            self.buffer = self.stream.read(self.bufsize)
            if not self.buffer: break
            self.pos = self.buffer.find('>')
            
    def __iter__(self):
        return self

    def next(self):
        if self.pos == -1:
            raise StopIteration
        
        startpos = self.pos + 1
        text = ""

        while 1:
            self.pos = self.buffer.find('>', startpos)

            if self.pos != -1:
                text += self.buffer[startpos: self.pos]
                break

            text += self.buffer[startpos:]
            self.buffer = self.stream.read(self.bufsize)

            if not self.buffer:
                break

            startpos = 0

        return self.parseRecord(text)

    def parseRecord(self, text):
        rec = Record()
        titlesplit = text.find('\n')
        rec.title = text[:titlesplit].rstrip()
        rec.sequence = text[titlesplit:].replace('\n', '').replace(' ', '')
        return rec

if __name__ == "__main__":
    import sys
    for rec in FastaIterator(file(sys.argv[1])):
        print rec
