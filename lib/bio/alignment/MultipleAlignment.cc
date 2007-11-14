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

#include <algorithm>

#include "bio/alignment/MultipleAlignment.hh"

namespace bio { namespace alignment {

	char MultipleAlignment::getChar(const size_t seqNum,
									const size_t column) const {
		return getSeq(seqNum)[column];
	}
	
	std::string MultipleAlignment::getSubstring(const size_t seqNum,
												const size_t colStart,
												const size_t colEnd) const {
		return getSeq(seqNum).substr(colStart, colEnd - colStart);
	}

	size_t MultipleAlignment::getSeqLen(const size_t seqNum) const {
		return getSeqLen(seqNum, 0, getNumCols());
	}
	
	size_t MultipleAlignment::getSeqLen(const size_t seqNum,
										const size_t colStart,
										const size_t colEnd) const {
		const std::string seq = getSeq(seqNum);
		size_t numGaps = std::count(seq.begin() + colStart,
									seq.begin() + colEnd,
									'-');
		return seq.size() - numGaps;
	}

	bool MultipleAlignment::isColIdentical(const size_t colNum) const {
		for (size_t i = 1; i < getNumSeqs(); ++i) {
			if (getChar(i, colNum) != getChar(0, colNum)) { return false; }
		}
		return true;
	}

	void MultipleAlignment::getGapMask(const size_t colNum,
									   std::vector<bool>& mask) const {
		mask.resize(getNumSeqs());
		for (size_t i = 0; i < getNumSeqs(); ++i) {
			mask[i] = (getChar(i, colNum) == '-');
		}
	}
	
	size_t MultipleAlignment::getNumBlocks() const {
		size_t num_blocks = 0;
		std::vector<bool> prev_mask, curr_mask;
		for (size_t col_num = 0; col_num < getNumCols(); ++col_num) {
			getGapMask(col_num, curr_mask);
			if (col_num == 0 or curr_mask != prev_mask) { ++num_blocks; }
			prev_mask.swap(curr_mask);
		}
		return num_blocks;
	}
	
	MultipleAlignment::~MultipleAlignment() {
	}
	
} }
