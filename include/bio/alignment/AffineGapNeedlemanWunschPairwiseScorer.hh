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

#ifndef __AFFINE_GAP_NEEDLEMAN_WUNSCH_PAIRWISE_SCORER_HH__
#define __AFFINE_GAP_NEEDLEMAN_WUNSCH_PAIRWISE_SCORER_HH__

#include <vector>

#include "bio/alignment/PairwiseSequenceScorer.hh"
#include "bio/alignment/ScoringMatrix.hh"

#include "util/stl.hh"
#include "util/matrix.hh"

namespace bio { namespace alignment {

	template<typename SemiRing>
	class AffineGapNeedlemanWunschPairwiseScorer
		: public PairwiseSequenceScorer<typename SemiRing::Element> {
	public:
		typedef typename SemiRing::Element Score;
		
		AffineGapNeedlemanWunschPairwiseScorer(const SemiRing& semiRing,
											   const ScoringMatrix<Score>& matrix,
											   const Score& space,
											   const Score& gap);
		
		Score score(const std::string& seq1,
					const std::string& seq2) const;

		static const size_t MATCH;
		static const size_t GAP1;
		static const size_t GAP2;
		
		typedef typename std::vector<Score> ScoreList;
		typedef typename util::Matrix<Score> ScoreMatrix;

		void scoreLastRow(const std::string& seq1,
						  const std::string& seq2,
						  ScoreList& matchRow,
						  ScoreList& gap1Row,
						  ScoreList& gap2Row,
						  size_t startState = MATCH) const;

		void scoreFirstRow(const std::string& seq1,
						   const std::string& seq2,
						   ScoreList& matchRow,
						   ScoreList& gap1Row,
						   ScoreList& gap2Row,
						   size_t endState = (MATCH | GAP1 | GAP2)) const;

		void scoreMatrixForward(const std::string& seq1,
								const std::string& seq2,
								ScoreMatrix& matchState,
								ScoreMatrix& gap1State,
								ScoreMatrix& gap2State) const;
		
		void scoreMatrixBackward(const std::string& seq1,
								 const std::string& seq2,
								 ScoreMatrix& matchState,
								 ScoreMatrix& gap1State,
								 ScoreMatrix& gap2State) const;

		void scoreMatrix(const std::string& seq1,
						 const std::string& seq2,
						 ScoreMatrix& matchState,
						 ScoreMatrix& gap1State,
						 ScoreMatrix& gap2State) const;
		
	protected:
		SemiRing semiRing;
		const ScoringMatrix<Score>& matrix;
		Score space;
		Score gap;
	};

	template<typename SemiRing>
	const size_t AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::MATCH = 1;

	template<typename SemiRing>
	const size_t AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::GAP1 = 2;

	template<typename SemiRing>
	const size_t AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::GAP2 = 4;

	template<typename SemiRing>
	AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::
	AffineGapNeedlemanWunschPairwiseScorer(const SemiRing& semiRing,
										   const ScoringMatrix<typename SemiRing::Element>& matrix,
										   const typename SemiRing::Element& space,
										   const typename SemiRing::Element& gap)
		: semiRing(semiRing), matrix(matrix), space(space), gap(gap) {
	}

	template<typename SemiRing>
	void
	AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::
	scoreLastRow(const std::string& seq1,
				 const std::string& seq2,
				 ScoreList& fromMatch,
				 ScoreList& fromGap1,
				 ScoreList& fromGap2,
				 size_t startState) const {
		const Score zero = semiRing.getZero();
		const Score one = semiRing.getMultiplicativeIdentity();

		// Initialize vectors
		fromMatch.resize(seq2.size() + 1, zero);
		fromGap1.resize(seq2.size() + 1, zero);
		fromGap2.resize(seq2.size() + 1, zero);

		fromMatch[0] = (startState & MATCH ? one : zero);
		fromGap1[0] = (startState & GAP1 ? one : zero);
		fromGap2[0] = (startState & GAP2 ? one : zero);

		// Calculate first row
		if (seq2.size() > 0) {
			fromGap2[1] = space * (gap * (fromMatch[0] + fromGap1[0]) +
								   fromGap2[0]);
			for (size_t j = 2; j <= seq2.size(); ++j) {
				fromGap2[j] = space * fromGap2[j - 1];
			}
		}

		// Calculate all rows
		Score upMatch, upGap1, upGap2, diagMatch, diagGap1, diagGap2;
		for (size_t i = 1; i <= seq1.size(); ++i) {
			for (size_t j = 0; j <= seq2.size(); ++j) {
				//std::cerr << "Node (" << i << "," << j << ")\n";
				
				upMatch = fromMatch[j];
				upGap1 = fromGap1[j];
				upGap2 = fromGap2[j];
				
				if (j == 0) {
					fromMatch[0] = fromGap2[0] = zero;
				} else {
					//std::cerr << "Match\n";
					fromMatch[j] = matrix.getCharScore(seq1[i - 1], seq2[j - 1]) *
						(diagMatch + diagGap1 + diagGap2);
					//std::cerr << "Gap2\n";
					fromGap2[j] = 
						space * (gap * (fromMatch[j - 1] + fromGap1[j - 1])
								 + fromGap2[j - 1]);
				}

				//std::cerr << "Gap1\n";
				fromGap1[j] = space * (gap * (upMatch + upGap2) + upGap1);
			
				diagMatch = upMatch;
				diagGap1 = upGap1;
				diagGap2 = upGap2;

			}
		}
	}
	
	template<typename SemiRing>
	void
	AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::
	scoreFirstRow(const std::string& seq1,
				  const std::string& seq2,
				  ScoreList& fromMatch,
				  ScoreList& fromGap1,
				  ScoreList& fromGap2,
				  size_t endState) const {
		const Score zero = semiRing.getZero();
		const Score one = semiRing.getMultiplicativeIdentity();
		
		// Initialize vectors
		fromMatch.resize(seq2.size() + 1, zero);
		fromGap1.resize(seq2.size() + 1, zero);
		fromGap2.resize(seq2.size() + 1, zero);

		fromMatch[seq2.size()] = (endState & MATCH ? one : zero);
		fromGap1[seq2.size()] = (endState & GAP1 ? one : zero);
		fromGap2[seq2.size()] = (endState & GAP2 ? one : zero);
		
		// Calculate last row
		for (size_t j = 0; j < seq2.size(); ++j) {
			size_t s = seq2.size() - j; // source
			size_t t = s - 1; // target
			fromGap2[t] = fromGap2[s] * space;
			fromMatch[t] = fromGap1[t] = fromGap2[t] * gap;
		}

		// Calculate all rows
		Score downMatch, diagMatch;
		for (size_t i = 1; i <= seq1.size(); ++i) {
			size_t tRow = seq1.size() - i;
			for (size_t j = 0; j <= seq2.size(); ++j) {
				size_t tCol = seq2.size() - j;
				size_t sCol = tCol + 1;
				
				downMatch = fromMatch[tCol];

				fromGap1[tCol] *= space;
				fromMatch[tCol] = fromGap2[tCol] = fromGap1[tCol] * gap;
				
				if (tCol < seq2.size()) {
					Score diag =
						diagMatch * matrix.getCharScore(seq1[tRow], seq2[tCol]);
					fromMatch[tCol] = fromGap2[tCol] += diag;
					fromGap1[tCol] += diag;

					fromMatch[tCol] += fromGap2[sCol] * gap * space;
					fromGap1[tCol] += fromGap2[sCol] * gap * space;
					fromGap2[tCol] += fromGap2[sCol] * space;
				}
			
				diagMatch = downMatch;
			}
		}
	}

	template<typename SemiRing>
	void
	AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::
	scoreMatrixForward(const std::string& seq1,
					   const std::string& seq2,
					   ScoreMatrix& matchState,
					   ScoreMatrix& gap1State,
					   ScoreMatrix& gap2State) const {
		const Score zero = semiRing.getZero();
		const Score one = semiRing.getMultiplicativeIdentity();

		// Initialize matrix
		matchState.resize(seq1.size() + 1, seq2.size() + 1);
		gap1State.resize(seq1.size() + 1, seq2.size() + 1);
		gap2State.resize(seq1.size() + 1, seq2.size() + 1);
		
		matchState(0, 0) = one;
		gap1State(0, 0) = zero;
		gap2State(0, 0) = zero;
		
		// Calculate first row
		if (seq2.size() > 0) {
			gap2State(0, 1) = gap * space;
			gap1State(0, 1) = matchState(0, 1) = zero;
			for (size_t j = 2; j <= seq2.size(); ++j) {
				gap2State(0, j) = space * gap2State(0, j - 1);
				gap1State(0, j) = matchState(0, j) = zero;
			}
		}

		// Calculate all rows
		for (size_t i = 1; i <= seq1.size(); ++i) {
			for (size_t j = 0; j <= seq2.size(); ++j) {
				if (j == 0) {
					matchState(i, j) = gap2State(i, j) = zero;
				} else {
					matchState(i, j) = (matchState(i - 1, j - 1) +
										gap1State(i - 1, j - 1) +
										gap2State(i - 1, j - 1)) *
						matrix.getCharScore(seq1[i - 1], seq2[j - 1]);
					gap2State(i, j) = space * (gap * (matchState(i, j - 1) +
													  gap1State(i, j - 1))
											   + gap2State(i, j - 1));
				}

				gap1State(i, j) = space * (gap * (matchState(i - 1, j) +
												  gap2State(i - 1, j)) +
										   gap1State(i - 1, j));
			}
		}
	}
		
	
	template<typename SemiRing>
	void
	AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::
	scoreMatrixBackward(const std::string& seq1,
						const std::string& seq2,
						ScoreMatrix& matchState,
						ScoreMatrix& gap1State,
						ScoreMatrix& gap2State) const {
		//const Score zero = semiRing.getZero();
		const Score one = semiRing.getMultiplicativeIdentity();

		size_t n = seq1.size();
		size_t m = seq2.size();
		
		// Initialize matrix
		matchState.resize(n + 1, m + 1);
		gap1State.resize(n + 1, m + 1);
		gap2State.resize(n + 1, m + 1);
		
		matchState(n, m) = gap1State(n, m) = gap2State(n, m) = one;

		// Calculate last row
		for (size_t j = 0; j < m; ++j) {
			size_t s = m - j; // source column
			size_t t = s - 1; // target column
			gap2State(n, t) = gap2State(n, s) * space;
			matchState(n, t) = gap1State(n, t) = gap2State(n, t) * gap;
		}

		// Calculate all rows
		for (size_t i = 1; i <= n; ++i) {
			size_t tRow = n - i; // target row
			size_t sRow = tRow + 1; // source row
			for (size_t j = 0; j <= m; ++j) {
				size_t tCol = m - j; // target column
				size_t sCol = tCol + 1; // source column

				gap1State(tRow, tCol) = gap1State(sRow, tCol) * space;
				matchState(tRow, tCol) = gap2State(tRow, tCol)
					= gap1State(tRow, tCol) * gap;
				
				if (tCol < m) {
					Score diag = matchState(sRow, sCol) * 
						matrix.getCharScore(seq1[tRow], seq2[tCol]);
					matchState(tRow, tCol) = gap2State(tRow, tCol) += diag;
					gap1State(tRow, tCol) += diag;

					Score horiz = gap2State(tRow, sCol) * space;
					matchState(tRow, tCol) += horiz * gap;
					gap1State(tRow, tCol) += horiz * gap;
					gap2State(tRow, tCol) += horiz;
				}
			}
		}
	}		

	template<typename SemiRing>
	void
	AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::
	scoreMatrix(const std::string& seq1,
				const std::string& seq2,
				ScoreMatrix& matchState,
				ScoreMatrix& gap1State,
				ScoreMatrix& gap2State) const {
		scoreMatrixForward(seq1, seq2, matchState, gap1State, gap2State);
		ScoreMatrix bMatchState, bGap1State, bGap2State;
		scoreMatrixBackward(seq1, seq2, bMatchState, bGap1State, bGap2State);
		matchState *= bMatchState;
		gap1State *= bGap1State;
		gap2State *= bGap2State;
	}
	
	template<typename SemiRing>
	typename SemiRing::Element
	AffineGapNeedlemanWunschPairwiseScorer<SemiRing>::
	score(const std::string& seq1,
		  const std::string& seq2) const {
		ScoreList fromMatch, fromGap1, fromGap2;
		scoreLastRow(seq1, seq2, fromMatch, fromGap1, fromGap2, MATCH);
		return fromMatch.back() + fromGap1.back() + fromGap2.back();
	}
		
} }

#endif // __AFFINE_GAP_NEEDLEMAN_WUNSCH_PAIRWISE_SCORER_HH__
