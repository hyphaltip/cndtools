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

#ifndef __AFFINE_GAP_SMITH_WATERMAN_PAIRWISE_SCORER_HH__
#define __AFFINE_GAP_SMITH_WATERMAN_PAIRWISE_SCORER_HH__

#include <vector>

#include "bio/alignment/PairwiseSequenceScorer.hh"

namespace bio { namespace alignment {

	template<typename SemiRing>
	class AffineGapSmithWatermanPairwiseScorer
		: public PairwiseSequenceScorer<typename SemiRing::Element> {
	public:
		typedef typename SemiRing::Element Score;
		
		AffineGapSmithWatermanPairwiseScorer(const Score& match,
											 const Score& mismatch,
											 const Score& space,
											 const Score& gap);
		
		Score score(const std::string& seq1,
					const std::string& seq2) const;

	protected:
		typedef typename std::vector<Score> ScoreList;

	private:
		Score match;
		Score mismatch;
		Score space;
		Score gap;
	};

	template<typename SemiRing>
	AffineGapSmithWatermanPairwiseScorer<SemiRing>::
	AffineGapSmithWatermanPairwiseScorer(const typename SemiRing::Element& match,
										 const typename SemiRing::Element& mismatch,
										 const typename SemiRing::Element& space,
										 const typename SemiRing::Element& gap)
		: match(match), mismatch(mismatch),
		  space(space), gap(gap) {
	}
		
	template<typename SemiRing>
	typename SemiRing::Element
	AffineGapSmithWatermanPairwiseScorer<SemiRing>::	
	score(const std::string& seq1,
		  const std::string& seq2) const {
		const Score zero = SemiRing::zero;
		const Score one = SemiRing::multiplicativeIdentity;

		// Initialize vectors
		std::vector<Score> fromMatch, fromGap1, fromGap2;
		fromMatch = fromGap1 = fromGap2 =
			std::vector<Score>(seq2.size() + 1, zero);
		fromMatch[0] = one;

		// Initialize result
		Score result = fromMatch[0];

		// Calculate first row
		if (seq2.size() > 0) {
			fromGap2[1] = gap * space;
			result += fromGap2[1];
			for (size_t j = 2; j <= seq2.size(); ++j) {
				fromGap2[j] = space * fromGap2[j - 1];
				result += fromGap2[j];
			}
		}
	
		// Calculate all rows
		Score upMatch, upGap1, upGap2, diagMatch, diagGap1, diagGap2;
		for (size_t i = 1; i <= seq1.size(); ++i) {
			for (size_t j = 0; j <= seq2.size(); ++j) {
				upMatch = fromMatch[j];
				upGap1 = fromGap1[j];
				upGap2 = fromGap2[j];

				if (j == 0) {
					fromMatch[0] = fromGap2[0] = zero;
				} else {
					fromMatch[j] =
						(seq1[i - 1] == seq2[j - 1] ? match : mismatch) *
						(diagMatch + diagGap1 + diagGap2);
					fromGap2[j] =
						space * (gap * (fromMatch[j - 1] + fromGap1[j - 1]) +
								 fromGap2[j - 1]);
					result += fromMatch[j];
					result += fromGap2[j];
				}

				fromGap1[j] = space * (gap * (upMatch + upGap2) + upGap1);
				result += fromGap1[j];
				
				diagMatch = upMatch;
				diagGap1 = upGap1;
				diagGap2 = upGap2;
			}
		}

		return result;
	}

	
} }

#endif // __AFFINE_GAP_SMITH_WATERMAN_PAIRWISE_SCORER_HH__
