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

#ifndef __AFFINE_GAP_PAIRWISE_ALIGNMENT_SCORER_HH__
#define __AFFINE_GAP_PAIRWISE_ALIGNMENT_SCORER_HH__

#include "bio/alignment/PairwiseAlignmentScorer.hh"
#include "bio/alignment/ScoringMatrix.hh"

namespace bio { namespace alignment {

	template<typename Score>
	class AffineGapPairwiseAlignmentScorer
		: public PairwiseAlignmentScorer<Score> {
	public:
		AffineGapPairwiseAlignmentScorer(const ScoringMatrix<Score>& matrix,
										 const Score& space,
										 const Score& gap);
		
		Score score(const PairwiseAlignment& pa) const;

	private:
		const ScoringMatrix<Score>& matrix;		
		Score space;
		Score gap;
	};

	template<typename Score>
	AffineGapPairwiseAlignmentScorer<Score>::
	AffineGapPairwiseAlignmentScorer(const ScoringMatrix<Score>& matrix,
									 const Score& space,
									 const Score& gap)
		: matrix(matrix), space(space), gap(gap) {
	}
	
	template<typename Score>
	Score
	AffineGapPairwiseAlignmentScorer<Score>::
    score(const PairwiseAlignment& pa) const {
		bool inGap1 = false;
		bool inGap2 = false;
		Score result = 0;
		for (std::string::size_type i = 0; i < pa.length(); ++i) {
			if (pa.seq1[i] == '-') {
				// Skip over positions that are gaps in both sequences
				// (useful for pairwise scoring of multiple
				// alignments)
				if (pa.seq2[i] == '-') {
					continue;
				}
				result += (inGap1 ? space : gap + space);
				inGap1 = true;
				inGap2 = false;
			} else if (pa.seq2[i] == '-') {
				result += (inGap2 ? space : gap + space);
				inGap1 = false;
				inGap2 = true;
			} else {
				result += matrix.getCharScore(pa.seq1[i], pa.seq2[i]);
				inGap1 = false;
				inGap2 = false;
			}
		}
		return result;
	}
	
} }

#endif // __AFFINE_GAP_PAIRWISE_ALIGNMENT_SCORER_HH__
