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
#include <cassert>

#include "bio/formats/phylip/OutputStream.hh"

namespace bio { namespace formats { namespace phylip {

	OutputStream::OutputStream(std::ostream& strm)
		: strm(strm),
		  hasWrittenRecord(false),
		  sequential(DEFAULT_SEQUENTIAL),
		  lineWidth(DEFAULT_LINE_WIDTH),
		  spaceFreq(DEFAULT_SPACE_FREQ),
		  blockSpacing(DEFAULT_BLOCK_SPACING)
	{}

	void OutputStream::setSequential(bool b) {
		this->sequential = b;
	}

	void OutputStream::setLineWidth(size_t lineWidth) {
		this->lineWidth = lineWidth;
	}

	void OutputStream::setSpaceFreq(size_t spaceFreq) {
		this->spaceFreq = spaceFreq;
	}

	void OutputStream::setBlockSpacing(size_t blockSpacing) {
		this->blockSpacing = blockSpacing;
	}

	void OutputStream::writeBlankLines(size_t numLines) {
		for (size_t i = 0; i < numLines; ++i) {
			strm << '\n';
		}
	}
	
	void OutputStream::writeFirstLine(const Record& rec) {
		strm << '\t' << rec.getNumSeqs()
			 << '\t' << rec.getNumCols() << '\n';
	}

	size_t OutputStream::writeSeq(const std::string& s,
								  size_t colsWritten,
								  size_t lineWidthRemaining) {
		assert(lineWidthRemaining > 0 and colsWritten < s.size());
		while (true) {

			size_t colsToWrite = std::min(s.size() - colsWritten,
										  lineWidthRemaining);
			if (spaceFreq > 0) {
				colsToWrite = std::min(colsToWrite,
									   spaceFreq - (colsWritten % spaceFreq));
			}
			strm << s.substr(colsWritten, colsToWrite);
			colsWritten += colsToWrite;
			lineWidthRemaining -= colsToWrite;
			if (colsWritten < s.size() and lineWidthRemaining > 0) {
				strm << ' ';
				--lineWidthRemaining;
			} else {
				break;
			}
		}
		strm << '\n';
		return colsWritten;
	}
	
	size_t OutputStream::writeSeqLine(const Record& rec,
									  size_t i,
									  size_t colsWritten) {
		//		strm << std::setw(MAX_NAME_LENGTH) << "";
		return writeSeq(rec.sequences[i], colsWritten, lineWidth);
		//						lineWidth - MAX_NAME_LENGTH);
	}
	
	size_t OutputStream::writeFirstSeqLine(const Record& rec, size_t i) {
		strm << std::setw(MAX_NAME_LENGTH)
			 << std::left
			 << rec.names[i].substr(0, MAX_NAME_LENGTH)
			 << ' ';
		return writeSeq(rec.sequences[i], 0,
						lineWidth - MAX_NAME_LENGTH);
	}
	
	void OutputStream::writeInterleaved(const Record& rec) {
		size_t colsWritten = 0;
		for (size_t i = 0; i < rec.getNumSeqs(); ++i) {
			colsWritten = writeFirstSeqLine(rec, i);
		}
		while (colsWritten < rec.getNumCols()) {
			writeBlankLines(blockSpacing);
			size_t nextColsWritten = 0;
			for (size_t i = 0; i < rec.getNumSeqs(); ++i) {
				nextColsWritten = writeSeqLine(rec, i, colsWritten);
			}
			colsWritten = nextColsWritten;
		}
	}

	void OutputStream::writeSequential(const Record& rec) {
		size_t colsWritten = 0;
		for (size_t i = 0; i < rec.getNumSeqs(); ++i) {
			colsWritten = writeFirstSeqLine(rec, i);
			while (colsWritten < rec.getNumCols()) {
				colsWritten = writeSeqLine(rec, i, colsWritten);
			}
		}
	}
	
	OutputStream& OutputStream::operator<<(const Record& rec) {
		if (hasWrittenRecord) {
			writeBlankLines(blockSpacing);
		} else {
			hasWrittenRecord = true;
		}
		writeFirstLine(rec);
		if (sequential) {
			writeSequential(rec);
		} else {
			writeInterleaved(rec);
		}
		return *this;
	}
	
	OutputStream::operator bool() const { return strm; }
	bool OutputStream::operator!() const { return not strm; }

} } }
