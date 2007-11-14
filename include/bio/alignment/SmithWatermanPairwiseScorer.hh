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

#ifndef __SMITH_WATERMAN_PAIRWISE_SCORER_HH__
#define __SMITH_WATERMAN_PAIRWISE_SCORER_HH__

#include <vector>

#include "bio/alignment/PairwiseSequenceScorer.hh"

namespace bio { namespace alignment {

	template<typename SemiRing>
	class SmithWatermanPairwiseScorer
		: public PairwiseSequenceScorer<typename SemiRing::Element> {
	public:
		typedef typename SemiRing::Element Element;
		
		SmithWatermanPairwiseScorer(const Element& match,
									const Element& mismatch,
									const Element& space);
		
		Element score(const std::string& seq1,
					  const std::string& seq2) const;

	private:
		typedef typename std::vector<Element> ElementList;
		Element match;
		Element mismatch;
		Element space;
	};

	template<typename SemiRing>
	SmithWatermanPairwiseScorer<SemiRing>::
	SmithWatermanPairwiseScorer(const typename SemiRing::Element& match,
								const typename SemiRing::Element& mismatch,
								const typename SemiRing::Element& space)
		: match(match), mismatch(mismatch), space(space) {
	}
		
	template<typename SemiRing>
	typename SemiRing::Element
	SmithWatermanPairwiseScorer<SemiRing>::	
	score(const std::string& seq1,
		  const std::string& seq2) const {
		const Element one = SemiRing::multiplicativeIdentity;
		
		// Initialize vectors
		ElementList currRow(seq2.size() + 1);
		currRow[0] = one;

		// Initialize result
		Element result = currRow[0];

		// Calculate first row
		for (size_t j = 1; j <= seq2.size(); ++j) {
			currRow[j] = space * currRow[j - 1] + one;
			result += currRow[j];
		}
	
		// Calculate all rows
		Element up, diag;
		for (size_t i = 1; i <= seq1.size(); ++i) {
			diag = currRow[0];
			currRow[0] *= space;
			currRow[0] += one;
			result += currRow[0];
			for (size_t j = 1; j <= seq2.size(); ++j) {
				up = currRow[j];

				if (seq1[i - 1] == seq2[j - 1]) {
					currRow[j] = match * diag + space * (up + currRow[j - 1]);
				} else {
					currRow[j] = mismatch * diag + space * (up + currRow[j - 1]);
				}
				currRow[j] += one;
				result += currRow[j];

				diag = up;
			}
		}
	
		return result;
	}
	
} }

#endif // __SMITH_WATERMAN_PAIRWISE_SCORER_HH__
