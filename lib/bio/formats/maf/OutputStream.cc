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
#include <limits>
#include <iomanip>

#include "bio/formats/maf/OutputStream.hh"

namespace bio { namespace formats { namespace maf {
	
	OutputStream::OutputStream(std::ostream& strm, const Header& header)
		: strm(strm) {
		writeHeader(header);
	}

	OutputStream& OutputStream::operator<<(const Record& rec) {
		// Determine the maximum length of the fields
		size_t maxSrcLength = 0;
		genome::Position maxStart = 0;
		size_t maxSize = 0;
		genome::Distance maxSrcSize = 0;
		for (size_t i = 0; i < rec.sequences.size(); ++i) {
			const Sequence& seq = rec.sequences[i];
			maxSrcLength = std::max(maxSrcLength, seq.src.size());
			maxStart = std::max(maxStart, seq.start);
			maxSize = std::max(maxSize, seq.size);
			maxSrcSize = std::max(maxSrcSize, seq.srcSize);
		}
		size_t maxStartLength = numDigits(maxStart);
		size_t maxSizeLength = numDigits(maxSize);
		size_t maxSrcSizeLength = numDigits(maxSrcSize);

		// Write blank line between alignments
		strm << '\n';

		// Write alignment line
		strm << ALIGNMENT_LINE_PREFIX;
		writeVariables(rec);
		strm << '\n';
		
		// Write sequence lines
		for (size_t i = 0; i < rec.sequences.size(); ++i) {
			writeSequence(rec.sequences[i],
						  maxSrcLength,
						  maxStartLength,
						  maxSizeLength,
						  maxSrcSizeLength);
		}
		
		return *this;
	}
	
	OutputStream& OutputStream::operator<<(const std::string& s) {
		strm << COMMENT_LINE_PREFIX << s << '\n';
		return *this;
	}
	
	void OutputStream::writeHeader(const Header& h) {
		strm << HEADER_LINE_PREFIX;
		writeVariables(h);
		strm << '\n';
		*this << h.parameters;
	}
	
	void OutputStream::writeVariables(const VariableMap& m) {
		VariableMap::Map::const_iterator var;
		for (var = m.variables.begin(); var != m.variables.end(); ++var) {
			strm << ' ' << var->first << '=' << var->second;
		}
	}

	void OutputStream::writeSequence(const Sequence& s,
									 size_t maxSrcLength,
									 size_t maxStartLength,
									 size_t maxSizeLength,
									 size_t maxSrcSizeLength) {
		strm << SEQ_LINE_PREFIX
			 << ' '
			 << std::setw(maxSrcLength)
			 << std::right				
			 << s.src
			 << ' '
			 << std::setw(maxStartLength)
			 << std::right
			 << s.start
			 << ' '
			 << std::setw(maxSizeLength)
			 << s.size
			 << ' '
			 << s.strand
			 << ' '
			 << std::setw(maxSrcSizeLength)
			 << std::right
			 << s.srcSize
			 << ' '
			 << s.text
			 << '\n';
	}

	OutputStream::operator bool() const { return strm; }
	bool OutputStream::operator!() const { return not strm; }
	
} } }
