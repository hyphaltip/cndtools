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

#ifndef __BIO_ALIGNMENT_SUMOFPAIRSSCORER_HH__
#define __BIO_ALIGNMENT_SUMOFPAIRSSCORER_HH__

#include "bio/alignment/MultipleAlignmentScorer.hh"

namespace bio { namespace alignment {

	template<typename Score>
	class SumOfPairsScorer : public MultipleAlignmentScorer<Score> {
	private:
		const PairwiseAlignmentScorer<Score>& pairwiseScorer;
		
	public:
		SumOfPairsScorer(const PairwiseAlignmentScorer<Score>& pairwiseScorer);
		
		Score score(const MultipleAlignment& pa) const;
	};

	template<typename Score>
	SumOfPairsScorer<Score>::
    SumOfPairsScorer(const PairwiseAlignmentScorer<Score>& pairwiseScorer)
		: pairwiseScorer(pairwiseScorer) {
	}
	
	template<typename Score>
	Score
	SumOfPairsScorer<Score>::
    score(const MultipleAlignment& ma) const {
		Score sum = 0;
		for (size_t i = 0; i < ma.getNumSeqs(); ++i) {
			for (size_t j = i + 1; j < ma.getNumSeqs(); ++j) {
				PairwiseAlignment pa(ma.getSeq(i), ma.getSeq(j));
				sum += pairwiseScorer.score(pa);
			}
		}
		return sum;
	}

} }

#endif // __BIO_ALIGNMENT_SUMOFPAIRSSCORER_HH__
