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

#include <sstream>
#include <istream>
#include <functional>

#include "bio/formats/phylip/InputStream.hh"
#include "util/string.hh"
#include "util/stl.hh"
using util::string::stripRight;
using util::string::IsWhitespace;

namespace bio { namespace formats { namespace phylip {

	struct EndOfFileException {};
	
	InputStream::InputStream(std::istream& strm)
		: strm(strm),
		  sequential(DEFAULT_SEQUENTIAL) {
	}

	void InputStream::setSequential(bool sequential) {
		this->sequential = sequential;
	}

	InputStream::operator bool() const { return strm; }

	bool InputStream::operator!() const { return not strm; }

	void InputStream::readUntilNextNonBlankLine() {
		while (strm >> line and line.empty());
		if (not strm) { throw EndOfFileException(); }
	}

	void InputStream::appendSequence(std::string::const_iterator start,
									 std::string::const_iterator end,
									 std::string& s) const {
		util::stl::copy_if(start, end, std::back_inserter(s),
						   std::not1(IsWhitespace()));
	}
	
	void InputStream::readFirstSeqLine(Record& rec, size_t i) {
		readUntilNextNonBlankLine();
		if (line.size() < MAX_NAME_LENGTH) {
			throw std::runtime_error("Initial sequence line does not "
									 "contain enough characters for "
									 "name:\n" + line);
		}
		rec.names[i] = stripRight(line.substr(0, MAX_NAME_LENGTH));
		appendSequence(line.begin() + MAX_NAME_LENGTH, line.end(),
					   rec.sequences[i]);
	}

	void InputStream::readSeqLine(Record& rec, size_t i) {
		readUntilNextNonBlankLine();
		appendSequence(line.begin(), line.end(), rec.sequences[i]);
	}

	void InputStream::readInterleaved(Record& rec, size_t numCols) {
		for (size_t i = 0; i < rec.sequences.size(); ++i) {
			readFirstSeqLine(rec, i);
		}
		while (rec.sequences.front().size() < numCols) {
			for (size_t i = 0; i < rec.sequences.size(); ++i) {
				readSeqLine(rec, i);
			}
		}
	}

	void InputStream::readSequential(Record& rec, size_t numCols) {
		for (size_t i = 0; i < rec.sequences.size(); ++i) {
			readFirstSeqLine(rec, i);
			while (rec.sequences[i].size() < numCols) {
				readSeqLine(rec, i);
			}
		}
	}

	void InputStream::readFirstLine(size_t& numSeqs, size_t& numCols) {
		readUntilNextNonBlankLine();
		std::istringstream firstLineStream(line);
		firstLineStream >> numSeqs >> numCols;
		if (not firstLineStream or numSeqs == 0 or numCols == 0) {
			throw std::runtime_error("Invalid first line:\n" + line);
		}
	}

	void InputStream::initRec(Record& rec, size_t numSeqs, size_t numCols) const {
		rec.names.resize(numSeqs);
		rec.sequences.resize(numSeqs);
		for (size_t i = 0; i < numSeqs; ++i) {
			rec.sequences[i].clear();
			rec.sequences[i].reserve(numCols);
		}
	}
	
	InputStream& InputStream::operator>>(Record& rec) {
		size_t numSeqs, numCols;
		try {
			readFirstLine(numSeqs, numCols);
		} catch (const EndOfFileException& e) {
			return *this;
		}
		
		initRec(rec, numSeqs, numCols);

		try {
			if (sequential) {
				readSequential(rec, numCols);
			} else {
				readInterleaved(rec, numCols);
			}
		} catch (const EndOfFileException& e) {
			throw std::runtime_error("Encountered end of file before end of "
									 "block");
		}

		return *this;
	}

} } }
