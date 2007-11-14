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

#include <cassert>

#include "bio/alignment/BasicMultipleAlignment.hh"
#include "util/string.hh"

namespace bio { namespace alignment {

	BasicMultipleAlignment::BasicMultipleAlignment() :
		seqs() {
	}
		
	size_t BasicMultipleAlignment::getNumCols() const {
		return seqs.front().size();
	}
			
	size_t BasicMultipleAlignment::getNumSeqs() const {
		return seqs.size();
	}

	std::string BasicMultipleAlignment::getSeq(const size_t seqNum) const {
		return seqs[seqNum];
	}

	char BasicMultipleAlignment::getChar(const size_t seqNum,
										 const size_t column) const {
		return seqs[seqNum][column];
	}
	
	std::string BasicMultipleAlignment::getSubstring(const size_t seqNum,
													 const size_t colStart,
													 const size_t colEnd) const {
		return seqs[seqNum].substr(colStart, colEnd - colStart);
	}


	void BasicMultipleAlignment::addSeq(const std::string& seq) {
		seqs.push_back(seq);
	}

	void BasicMultipleAlignment::clear() {
		seqs.clear();
	}
	
	size_t BasicMultipleAlignment::getColumnNum(const size_t seqNum,
												const size_t seqPos) const {
		assert(seqNum < getNumSeqs());
		size_t nextPos = 0;
		for (size_t col = 0; col < getNumCols(); ++col) {
			if (seqs[seqNum][col] != '-') {
				if (nextPos == seqPos) {
					return col;
				} else {
					++nextPos;
				}
			}
		}
		assert(false);
		return 0;
	}

	size_t BasicMultipleAlignment::getSeqPos(const size_t seqNum,
										const size_t col) const {
		assert(seqNum < getNumSeqs());
		size_t seqPos = 0;
		for (std::string::const_iterator pos = seqs[seqNum].begin();
			 pos != seqs[seqNum].begin() + col;
			 ++pos) {
			if (*pos != '-') {
				++seqPos;
			}
		}
		return seqPos;
	}
	
	Interval BasicMultipleAlignment::getColumnInterval(const size_t seqNum,
													   const size_t seqStart,
													   const size_t seqEnd) const {
		assert(seqNum < getNumSeqs());
		assert(seqEnd >= seqStart);
		Interval colInt;
		size_t nextPos = 0;
		for (size_t col = 0; col < getNumCols(); ++col) {
			if (seqs[seqNum][col] != '-') {
				if (nextPos == seqStart) {
					colInt.start = col;
				}
				if (nextPos == (seqEnd - 1)) {
					colInt.end = col + 1;
					return colInt;
				}
				++nextPos;
			}
		}
		assert(false);
		return colInt;
	}

	Interval BasicMultipleAlignment::getSeqInterval(const size_t seqNum,
											   const size_t colStart,
											   const size_t colEnd) const {
		assert(seqNum < getNumSeqs());
		assert(colEnd >= colStart);

		Interval seqInt;
		seqInt.start = getSeqPos(seqNum, colStart);
		seqInt.end = seqInt.start;

		for (std::string::const_iterator pos = seqs[seqNum].begin() + colStart;
			 pos != seqs[seqNum].begin() + colEnd;
			 ++pos) {
			if (*pos != '-') {
				++seqInt.end;
			}
		}

		return seqInt;
	}
	

	std::string BasicMultipleAlignment::getCharsAfter(const size_t seqNum,
												 const size_t colStart,
												 const size_t numChars) const {
		std::string first(numChars, '?');
		std::string::iterator firstPos = first.begin();
		for (std::string::const_iterator seqPos = seqs[seqNum].begin() + colStart;
			 seqPos != seqs[seqNum].end() && firstPos != first.end();
			 ++seqPos) {
			if (*seqPos != '-') {
				(*firstPos) = *seqPos;
				++firstPos;
			}
		}
		return first;
	}
	
	std::string BasicMultipleAlignment::getCharsBefore(const size_t seqNum,
												  const size_t colEnd,
												  const size_t numChars) const {
		std::string last(numChars, '?');
		std::string::reverse_iterator lastPos = last.rbegin();
		for (std::string::const_reverse_iterator seqPos =
				 seqs[seqNum].rbegin() + (seqs[seqNum].size() - colEnd);
			 seqPos != seqs[seqNum].rend() && lastPos != last.rend();
			 ++seqPos) {
			if (*seqPos != '-') {
				(*lastPos) = *seqPos;
				++lastPos;
			}
		}
		return last;
	}

	size_t BasicMultipleAlignment::getNumChars(const size_t seqNum,
											   const size_t colStart,
											   const size_t colEnd) const {
		size_t numChars = 0;
		for (std::string::const_iterator
				 seqPos = seqs[seqNum].begin() + colStart;
			 seqPos != seqs[seqNum].begin() + colEnd;
			 ++seqPos) {
			if (*seqPos != '-') {
				++numChars;
			}
		}
		return numChars;
	}

	void BasicMultipleAlignment::makeSeqsUppercase() {
		for (size_t i = 0; i < seqs.size(); ++i) {
			util::string::upperString(seqs[i]);
		}
	}

} }
