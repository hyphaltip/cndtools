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

#include "util/stl.hh"
#include "bio/alignment/PairwiseAlignment.hh"

namespace bio { namespace alignment {

	size_t PairwiseAlignment::getNumIdentities() const {
		return util::stl::matches(seq1.begin(), seq1.end(), seq2.begin());
	}
	
	PairwiseAlignment
	PairwiseAlignment::slice(size_t seqNum, size_t start, size_t end) const {
		std::string seq = (seqNum == 0 ? seq1 : seq2);
		size_t startCol = seq.find_first_not_of('-');
		while (start > 0) {
			startCol = seq.find_first_not_of('-', startCol + 1);
			--start;
			--end;
		}
		size_t endCol = startCol;
		while (end > 0) {
			if (end == 1) {
				++endCol;
			} else {
				endCol = seq.find_first_not_of('-', endCol + 1);
			}
			--end;
		}
		return PairwiseAlignment(seq1.substr(startCol, endCol - startCol),
								 seq2.substr(startCol, endCol - startCol));
	}
	
} }
