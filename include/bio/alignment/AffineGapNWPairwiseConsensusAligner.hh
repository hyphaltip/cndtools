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

#ifndef __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISECONSENSUSALIGNER_HH__
#define __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISECONSENSUSALIGNER_HH__

#include <vector>

#include "bio/alignment/PairwiseAligner.hh"
#include "bio/alignment/AffineGapNeedlemanWunschPairwiseScorer.hh"
#include "bio/alignment/ConvertingScoringMatrix.hh"
#include "math/MaxPlus.hh"
#include "util/string.hh"

namespace bio { namespace alignment {

	template<typename NumType>
	class AffineGapNWPairwiseConsensusAligner
		: public PairwiseAligner,
		  protected AffineGapNeedlemanWunschPairwiseScorer< math::MaxPlus<NumType> > {
	public:
		
		AffineGapNWPairwiseConsensusAligner(const ScoringMatrix<NumType>& matrix,
											const NumType& space,
											const NumType& gap);
		
		PairwiseAlignment align(const std::string& seq1,
								const std::string& seq2) const;
				
	protected:
		typedef typename AffineGapNeedlemanWunschPairwiseScorer< math::MaxPlus<NumType> >::Score Score;
		typedef typename AffineGapNeedlemanWunschPairwiseScorer< math::MaxPlus<NumType> >::ScoreMatrix ScoreMatrix;

		void findConsensusEdges(const ScoreMatrix& matchState,
								const ScoreMatrix& gap1State,
								std::vector<size_t>& alignment,
								std::vector<bool>& isAligned) const;

		PairwiseAlignment constructAlignment(const std::string& seq1,
											 const std::string& seq2,
											 const std::vector<size_t>& alignment1,
											 const std::vector<bool>& isAligned1,
											 const std::vector<size_t>& alignment2,
											 const std::vector<bool>& isAligned2) const;

		ConvertingScoringMatrix<NumType, Score> convertingMatrix;
	};

	template<typename NumType>
	AffineGapNWPairwiseConsensusAligner<NumType>::
	AffineGapNWPairwiseConsensusAligner(const ScoringMatrix<NumType>& matrix,
										const NumType& space,
										const NumType& gap)
		: AffineGapNeedlemanWunschPairwiseScorer< math::MaxPlus<NumType> >(math::MaxPlus<NumType>(),
																		   convertingMatrix,
																		   space,
																		   gap),
		  convertingMatrix(matrix) {
	}

	template<typename NumType>
	void
	AffineGapNWPairwiseConsensusAligner<NumType>::
	findConsensusEdges(const ScoreMatrix& matchState,
					   const ScoreMatrix& gap1State,
					   std::vector<size_t>& alignment,
					   std::vector<bool>& isAligned) const {
		size_t n = matchState.getNumRows() - 1;
		size_t m = matchState.getNumCols() - 1;
		
		// Initialize alignment vectors
		alignment.resize(n);
		isAligned.resize(n);		
		std::fill(alignment.begin(), alignment.end(), m);
		std::fill(isAligned.begin(), isAligned.end(), false);
		
		Score optimal = matchState(0, 0);

		// Determine if a gap path is optimal
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j <= m; ++j) {
				bool isGapPathOpt = (gap1State(i + 1, j) == optimal);
				bool isMatchPathOpt = (j < m and
									   matchState(i + 1, j + 1) == optimal);
				
				if ((isGapPathOpt and isMatchPathOpt) or
					(isMatchPathOpt and isAligned[i]) or
					(isGapPathOpt and alignment[i] != m)) {
					isAligned[i] = false;
					break;
				} else if (isGapPathOpt or isMatchPathOpt) {
					isAligned[i] = true;
					if (isMatchPathOpt) {
						alignment[i] = j;
					}
				}
			}
		}
	}

	template<typename NumType>
	PairwiseAlignment
	AffineGapNWPairwiseConsensusAligner<NumType>::
	constructAlignment(const std::string& seq1,
					   const std::string& seq2,
					   const std::vector<size_t>& alignment1,
					   const std::vector<bool>& isAligned1,
					   const std::vector<size_t>& alignment2,
					   const std::vector<bool>& isAligned2) const {		
		std::string alignString1, alignString2;
		size_t i = 0, j = 0;
		while (i < seq1.size() or j < seq2.size()) {
			if (i != seq1.size() and
				(j == seq2.size()
				 or (not isAligned1[i])
				 or alignment1[i] == seq2.size())) {
				alignString1.push_back(seq1.at(i));
				alignString2.push_back(isAligned1[i] ? '-' : '=');
				++i;
			} else if (i == seq1.size() or (not isAligned2[j])
					   or alignment2[j] == seq1.size()) {
				alignString1.push_back(isAligned2[j] ? '-' : '=');
				alignString2.push_back(seq2.at(j));
				++j;
			} else {
				alignString1.push_back(seq1.at(i));
				alignString2.push_back(seq2.at(j));
				++i;
				++j;
			}
		}
		
		return PairwiseAlignment(alignString1, alignString2);
	}
	
	template<typename NumType>
	PairwiseAlignment
	AffineGapNWPairwiseConsensusAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2) const {
		ScoreMatrix matchState, gap1State, gap2State;
		scoreMatrix(seq1, seq2, matchState, gap1State, gap2State);
		
		std::vector<bool> isAligned1, isAligned2;
		std::vector<size_t> alignment1, alignment2;

		findConsensusEdges(matchState, gap1State,
						   alignment1, isAligned1);
		matchState.transpose();
		gap1State.transpose();
		gap2State.transpose();
		findConsensusEdges(matchState, gap2State,
						   alignment2, isAligned2);
		
		return constructAlignment(seq1, seq2,
								  alignment1, isAligned1,
								  alignment2, isAligned2);
	}

} }

#endif // __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISECONSENSUSALIGNER_HH__
