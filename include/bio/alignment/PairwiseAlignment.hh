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

#ifndef __PAIRWISE_ALIGNMENT_HH__
#define __PAIRWISE_ALIGNMENT_HH__

#include <string>

namespace bio { namespace alignment {

	class PairwiseAlignment {
	public:
		std::string seq1;
		std::string seq2;

		PairwiseAlignment(const std::string& seq1 = "",
						  const std::string& seq2 = "")
			: seq1(seq1), seq2(seq2) {}

		void flip() { std::swap(seq1, seq2); }

		std::string::size_type length() const { return seq1.length(); }

		size_t getNumIdentities() const;

		PairwiseAlignment slice(size_t seqNum, size_t start, size_t end) const;
	   
		PairwiseAlignment& operator+=(const PairwiseAlignment& other) {
			seq1 += other.seq1;
			seq2 += other.seq2;
			return *this;
		}

		PairwiseAlignment operator+(const PairwiseAlignment& other) const {
			PairwiseAlignment pa = *this;
			return pa += other;
		}

		bool operator<(const PairwiseAlignment& other) const {
			return seq1 < other.seq1 or
				(seq1 == other.seq1 and seq2 < other.seq2);
		}

		bool operator==(const PairwiseAlignment& other) const {
			return seq1 == other.seq1 and seq2 == other.seq2;
		}

	};

} }

#endif // __PAIRWISE_ALIGNMENT_HH__
