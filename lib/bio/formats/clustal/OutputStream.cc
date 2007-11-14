/* Copyright (c) 2006
   Colin Dewey (University of Wisconsin-Madison)
   cdewey@biostat.wisc.edu
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <ostream>
#include <iomanip>
#include <stdexcept>

#include "bio/formats/clustal/OutputStream.hh"
#include "util/string.hh"
using util::string::firstWord;

namespace bio { namespace formats { namespace clustal {

	OutputStream::OutputStream(std::ostream& strm)
		: strm(strm),
		  lineWidth(DEFAULT_LINE_WIDTH),
		  minNameSeqSpacing(DEFAULT_MIN_NAME_SEQ_SPACING),
		  minGutterLen(DEFAULT_MIN_GUTTER_LEN),
		  blockSpacing(DEFAULT_BLOCK_SPACING),
		  seqnos(DEFAULT_SEQNOS)
	{}

	void OutputStream::setLineWidth(size_t lineWidth) {
		this->lineWidth = lineWidth;
	}

	void OutputStream::setMinNameSeqSpacing(size_t minNameSeqSpacing) {
		this->minNameSeqSpacing = minNameSeqSpacing;
	}

	void OutputStream::setMinGutterLen(size_t minGutterLen) {
		this->minGutterLen = minGutterLen;
	}

	void OutputStream::setBlockSpacing(size_t blockSpacing) {
		this->blockSpacing = blockSpacing;
	}

	void OutputStream::setSeqNos(bool b) {
		this->seqnos = b;
	}
	
	void OutputStream::writeBlankLines(size_t numLines) {
		for (size_t i = 0; i < numLines; ++i) {
			strm << '\n';
		}
	}
	
	void OutputStream::writeFirstLine() {
		strm << HEADER_LINE_PREFIX << '\n';
	}

	void OutputStream::writeBlockLine(const std::string& name,
									  const std::string& seq,
									  size_t gutterLen,
									  size_t& seqPos) {
		size_t numGaps = std::count(seq.begin(), seq.end(), '-');
		seqPos += seq.length() - numGaps;
		strm << std::setw(gutterLen) << std::left << name << seq;
		if (seqnos and seq.length() != numGaps) { strm << ' ' << seqPos; }
		strm << '\n';
	}

	void OutputStream::writeBlock(const Record& rec,
								  size_t start,
								  size_t gutterLen,
								  size_t seqLen,
								  std::vector<size_t>& seqPos) {
		writeBlankLines(1);
		for (size_t i = 0; i < rec.names.size(); ++i) {
			writeBlockLine(rec.names[i],
						   rec.sequences[i].substr(start, seqLen),
						   gutterLen,
						   seqPos[i]);
		}
		writeIdentityLine(rec.sequences, start, seqLen, gutterLen);
	}

	void OutputStream::writeIdentityLine(const std::vector<std::string>& v,
										 size_t start,
										 size_t seqLen,
										 size_t gutterLen) {
		strm << std::setw(gutterLen) << "";
		size_t stop = std::min(start + seqLen, v.front().size());
		for (size_t pos = start; pos < stop; ++pos) {
			strm << (isIdenticalColumn(v, pos) ? '*' : ' ');
		}
		strm << '\n';
	}

	bool OutputStream::isIdenticalColumn(const std::vector<std::string>& v,
										 size_t pos) {
		char c = v.front()[pos];
		for (size_t i = 1; i < v.size(); ++i) {
			if (v[i][pos] != c) { return false; }
		}
		return true;
	}

	size_t OutputStream::maxLen(const std::vector<std::string>& v) {
		size_t len = 0;
		for (size_t i = 0; i < v.size(); ++i) {
			len = std::max(len, v[i].size());
		}
		return len;
	}

	size_t OutputStream::numDigits(size_t n) {
		size_t digits = 1;
		while (n >= 10) {
			n /= 10;
			++digits;
		}
		return digits;
	}
	
	size_t OutputStream::maxSeqNosLen(const std::vector<std::string>& v) {
		size_t len = 0;
		for (size_t i = 0; i < v.size(); ++i) {
			len = std::max(len, numDigits(v[i].size()));
		}
		return len;
	}
	
	void OutputStream::writeBlocks(const Record& rec) {
		if (rec.sequences.empty()) { return; }
		size_t gutterLen = std::max(minGutterLen,
									maxLen(rec.names) + minNameSeqSpacing);
		if (lineWidth < gutterLen) {
			throw std::runtime_error("Gutter width is less than line width");
		}
		size_t seqNosLen = (seqnos ? 1 + maxSeqNosLen(rec.sequences) : 0);
		size_t seqLen = lineWidth - gutterLen - seqNosLen;
		size_t alignmentLen = rec.sequences.front().size();
		std::vector<size_t> seqPos(rec.getNumSeqs(), 0);
		for (size_t start = 0; start < alignmentLen; start += seqLen) {
			writeBlock(rec, start, gutterLen, seqLen, seqPos);
		}
	}
	
	OutputStream& OutputStream::operator<<(const Record& rec) {
		writeFirstLine();
		writeBlocks(rec);
		return *this;
	}

	OutputStream& OutputStream::
	operator<<(const alignment::NamedMultipleAlignment& align) {
		Record rec;
		for (size_t i = 0; i < align.getNumSeqs(); ++i) {
			rec.names.push_back(firstWord(align.getName(i)));
			rec.sequences.push_back(align.getSeq(i));
		}
		return *this << rec;
	}
	
	OutputStream::operator bool() const { return strm; }
	bool OutputStream::operator!() const { return not strm; }

} } }
