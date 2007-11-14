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

#ifndef __BIO_ALIGNMENT_AFFINEGAPNWFASTPAIRWISEALIGNER3_HH__
#define __BIO_ALIGNMENT_AFFINEGAPNWFASTPAIRWISEALIGNER3_HH__

#include "bio/alignment/PairwiseAligner.hh"
#include "bio/alignment/AffineGapNWPairwiseScorer.hh"
#include "bio/alignment/ConvertingScoringMatrix.hh"
#include "math/MaxPlus.hh"
#include "util/stl.hh"

namespace bio { namespace alignment {

	template<typename NumType>	
    class AffineGapNWFastPairwiseAligner3 : public PairwiseAligner {
	public:
		AffineGapNWFastPairwiseAligner3(const ScoringMatrix<NumType>& matrix,
										const NumType& space,
										const NumType& gap);
		
		PairwiseAlignment align(const std::string& seq1,
								const std::string& seq2) const;

	private:
		typedef typename math::MaxPlus<NumType> SemiRing;
		typedef AffineGapNWPairwiseScorer<SemiRing> Scorer;
  		typedef typename Scorer::Score Score;
  		typedef typename Scorer::ScoreMatrix ScoreMatrix;
		
		static PairwiseAlignment
		statesToAlignment(const std::vector<int>& indices,
						  const std::string& seq1,
						  const std::string& seq2);

		static int argmax(Score h, Score d, Score i);

		static const int H;
		static const int I;
		static const int D;

		ConvertingScoringMatrix<NumType, Score> convertingMatrix;
		Score space;
		Score gap;
		Scorer scorer;
	};

	template<typename NumType>	
    const int AffineGapNWFastPairwiseAligner3<NumType>::H = 0;
	template<typename NumType>	
    const int AffineGapNWFastPairwiseAligner3<NumType>::D = 1;
	template<typename NumType>	
    const int AffineGapNWFastPairwiseAligner3<NumType>::I = 2;
	
	template<typename NumType>
	AffineGapNWFastPairwiseAligner3<NumType>::
	AffineGapNWFastPairwiseAligner3(const ScoringMatrix<NumType>& matrix,
									const NumType& space,
									const NumType& gap)
		: convertingMatrix(matrix),
		  space(space),
		  gap(gap),
		  scorer(SemiRing(), convertingMatrix, space, gap) {
	}
	
	template<typename NumType>
	int
	AffineGapNWFastPairwiseAligner3<NumType>::
	argmax(Score h, Score d, Score i) {
		if (h < d and i < d) {
			return D;
		} else if (h < i) {
			return I;
		} else {
			return H;
		}
	}

	template<typename NumType>
	PairwiseAlignment
	AffineGapNWFastPairwiseAligner3<NumType>::
	statesToAlignment(const std::vector<int>& indices,
					  const std::string& seq1,
					  const std::string& seq2) {
		PairwiseAlignment alignment;
		size_t pos1 = 0, pos2 = 0;
		for (size_t i = 0; i < indices.size(); ++i) {
			if (indices[i] == H) {
				alignment.seq1 += seq1[pos1];
				alignment.seq2 += seq2[pos2];
				++pos1;
				++pos2;
			} else if (indices[i] == I) {
				alignment.seq1 += '-';
				alignment.seq2 += seq2[pos2];
				++pos2;
			} else {
				alignment.seq1 += seq1[pos1];
				alignment.seq2 += '-';
				++pos1;
			}
		}

		return alignment;
	}
	
	template<typename NumType>
	PairwiseAlignment
	AffineGapNWFastPairwiseAligner3<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2) const {
		std::vector<ScoreMatrix> states(3);
		scorer.scoreMatrixForward(seq1, seq1, states[H], states[D], states[I]);

		std::vector<int> traceback;

		size_t i = seq1.size(), j = seq2.size();
		
		// Figure out the best state to end in
		int currState = argmax(states[H](i, j),
							   states[D](i, j),
							   states[I](i, j));
		traceback.push_back(currState);

		while (i > 0 or j > 0) {
			if (currState == H) {
				--i;
				--j;
				currState = argmax(states[H](i, j),
								   states[D](i, j),
								   states[I](i, j));
			} else if (currState == D) {
				--i;
				currState = argmax(states[H](i, j) * gap,
								   states[D](i, j),
								   states[I](i, j) * gap);
			} else if (currState == I) {
				--j;
				currState = argmax(states[H](i, j) * gap,
								   states[D](i, j) * gap,
								   states[I](i, j));
			} else {
				throw std::runtime_error("AffineGapNWFastPairwiseAligner3"
										 "::align: Bad state");
			}
			traceback.push_back(currState);
		}
		std::reverse(traceback.begin(), traceback.end());
		
		return statesToAlignment(traceback, seq1, seq2);
	}
	
} }

#endif // __BIO_ALIGNMENT_AFFINEGAPNWFASTPAIRWISEALIGNER3_HH__
