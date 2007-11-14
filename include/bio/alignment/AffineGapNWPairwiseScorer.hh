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

#ifndef __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISESCORER_HH__
#define __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISESCORER_HH__

#include <vector>

#include "bio/alignment/PairwiseSequenceScorer.hh"
#include "bio/alignment/AlphabetScoringMatrix.hh"

#include "util/stl.hh"
#include "util/matrix.hh"

using util::stl::print_elements;

namespace bio { namespace alignment {

	template<typename SemiRing>
	class AffineGapNWPairwiseScorer
		: public PairwiseSequenceScorer<typename SemiRing::Element> {
	public:
		typedef typename SemiRing::Element Score;
		
		AffineGapNWPairwiseScorer(const SemiRing& semiRing,
								  const ScoringMatrix<Score>& matrix,
								  const Score& space,
								  const Score& gap);
		
		AffineGapNWPairwiseScorer(const SemiRing& semiRing,
								  const ScoringMatrix<Score>& matrix,
								  const Score& insertionSpace,
								  const Score& deletionSpace,
								  const Score& insertionGap,
								  const Score& deletionGap);

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
		AlphabetScoringMatrix<Score> matrix;
		Score insertionSpace;
		Score deletionSpace;
		Score insertionGap;
		Score deletionGap;
	};

	template<typename SemiRing>
	const size_t AffineGapNWPairwiseScorer<SemiRing>::MATCH = 1;

	template<typename SemiRing>
	const size_t AffineGapNWPairwiseScorer<SemiRing>::GAP1 = 2;

	template<typename SemiRing>
	const size_t AffineGapNWPairwiseScorer<SemiRing>::GAP2 = 4;

	template<typename SemiRing>
	AffineGapNWPairwiseScorer<SemiRing>::
	AffineGapNWPairwiseScorer(const SemiRing& semiRing,
							  const ScoringMatrix<Score>& matrix,
							  const Score& space,
							  const Score& gap)
		: semiRing(semiRing), matrix(matrix),
		  insertionSpace(space), deletionSpace(space),
		  insertionGap(gap), deletionGap(gap) {
	}

	template<typename SemiRing>
	AffineGapNWPairwiseScorer<SemiRing>::
	AffineGapNWPairwiseScorer(const SemiRing& semiRing,
							  const ScoringMatrix<Score>& matrix,
							  const Score& insertionSpace,
							  const Score& deletionSpace,
							  const Score& insertionGap,
							  const Score& deletionGap)
		: semiRing(semiRing), matrix(matrix),
		  insertionSpace(insertionSpace),
		  deletionSpace(deletionSpace),
		  insertionGap(insertionGap),
		  deletionGap(deletionGap) {
	}
	
	template<typename SemiRing>
	void
	AffineGapNWPairwiseScorer<SemiRing>::
	scoreLastRow(const std::string& seq1,
				 const std::string& seq2,
				 ScoreList& fromMatch,
				 ScoreList& fromDeletion,
				 ScoreList& fromInsertion,
				 size_t startState) const {
		const Score zero = semiRing.getZero();
		const Score one = semiRing.getMultiplicativeIdentity();

		size_t rows = seq1.size() + 1;
		size_t cols = seq2.size() + 1;
		
		// Initialize vectors
		fromMatch.resize(cols, zero);
		fromDeletion.resize(cols, zero);
		fromInsertion.resize(cols, zero);

		ScoreList upMatch(cols, zero);
		ScoreList upDeletion(cols, zero);
		ScoreList upInsertion(cols, zero);

		// Calculate first row
		fromMatch[0] = (startState & MATCH ? one : zero);
		fromInsertion[0] = (startState & GAP2 ? one : zero);
		fromDeletion[0] = (startState & GAP1 ? one : zero);		

		if (cols > 1) {
			fromInsertion[1] = fromMatch[0];
			fromInsertion[1] += fromDeletion[0];
			fromInsertion[1] *= insertionGap;
			fromInsertion[1] += fromInsertion[0];
			fromInsertion[1] *= insertionSpace;
			
			for (size_t j = 2; j < cols; ++j) {
				fromInsertion[j] = fromInsertion[j - 1];
				fromInsertion[j] *= insertionSpace;
			}
		}
		
		// Calculate all rows
		for (size_t i = 1; i < rows; ++i) {
			fromMatch.swap(upMatch);
			fromDeletion.swap(upDeletion);
			fromInsertion.swap(upInsertion);
			
			for (size_t j = 0; j < cols; ++j) {

				//std::cerr << "Node (" << i << "," << j << ")\n";
				
				if (j == 0) {
					fromMatch[0] = zero;
					fromInsertion[0] = zero;
				} else {
					//std::cerr << "Match\n";
					fromMatch[j] = upMatch[j - 1];
					fromMatch[j] += upDeletion[j - 1];
					fromMatch[j] += upInsertion[j - 1];
					fromMatch[j] *= matrix.getCharScore(seq1[i - 1], seq2[j - 1]);

					//std::cerr << "Insertion\n";
					fromInsertion[j] = fromMatch[j - 1];
					fromInsertion[j] += fromDeletion[j - 1];
					fromInsertion[j] *= insertionGap;
					fromInsertion[j] += fromInsertion[j - 1];
					fromInsertion[j] *= insertionSpace;

				}
				
				//std::cerr << "Deletion\n";
				fromDeletion[j] = upMatch[j];
				fromDeletion[j] += upInsertion[j];
				fromDeletion[j] *= deletionGap;
				fromDeletion[j] += upDeletion[j];
				fromDeletion[j] *= deletionSpace;
			}
		}
	}
	
	template<typename SemiRing>
	void
	AffineGapNWPairwiseScorer<SemiRing>::
	scoreFirstRow(const std::string& seq1,
				  const std::string& seq2,
				  ScoreList& fromMatch,
				  ScoreList& fromDeletion,
				  ScoreList& fromInsertion,
				  size_t endState) const {
		const Score zero = semiRing.getZero();
		const Score one = semiRing.getMultiplicativeIdentity();
		
		// Initialize vectors
		fromMatch.resize(seq2.size() + 1, zero);
		fromDeletion.resize(seq2.size() + 1, zero);
		fromInsertion.resize(seq2.size() + 1, zero);

		fromMatch[seq2.size()] = (endState & MATCH ? one : zero);
		fromDeletion[seq2.size()] = (endState & GAP1 ? one : zero);
		fromInsertion[seq2.size()] = (endState & GAP2 ? one : zero);
		
		// Calculate last row
		for (size_t j = 0; j < seq2.size(); ++j) {
			size_t s = seq2.size() - j; // source
			size_t t = s - 1; // target
			fromInsertion[t] = fromInsertion[s] * insertionSpace;
			fromMatch[t] = fromDeletion[t] = fromInsertion[t] * insertionGap;
		}

		// Calculate all rows
		Score downMatch, diagMatch;
		for (size_t i = 1; i <= seq1.size(); ++i) {
			size_t tRow = seq1.size() - i;
			for (size_t j = 0; j <= seq2.size(); ++j) {
				size_t tCol = seq2.size() - j;
				size_t sCol = tCol + 1;
				
				downMatch = fromMatch[tCol];

				fromDeletion[tCol] *= deletionSpace;
				fromMatch[tCol] = fromInsertion[tCol] = fromDeletion[tCol] * deletionGap;
				
				if (tCol < seq2.size()) {
					Score diag =
						diagMatch * matrix.getCharScore(seq1[tRow], seq2[tCol]);
					fromMatch[tCol] = fromInsertion[tCol] += diag;
					fromDeletion[tCol] += diag;

					fromMatch[tCol] += fromInsertion[sCol] * insertionGap * insertionSpace;
					fromDeletion[tCol] += fromInsertion[sCol] * insertionGap * insertionSpace;
					fromInsertion[tCol] += fromInsertion[sCol] * insertionSpace;
				}
			
				diagMatch = downMatch;
			}
		}
	}

	template<typename SemiRing>
	void
	AffineGapNWPairwiseScorer<SemiRing>::
	scoreMatrixForward(const std::string& seq1,
					   const std::string& seq2,
					   ScoreMatrix& matchState,
					   ScoreMatrix& deletionState,
					   ScoreMatrix& insertionState) const {
		const Score zero = semiRing.getZero();
		const Score one = semiRing.getMultiplicativeIdentity();

		// Initialize matrix
		matchState.resize(seq1.size() + 1, seq2.size() + 1);
		deletionState.resize(seq1.size() + 1, seq2.size() + 1);
		insertionState.resize(seq1.size() + 1, seq2.size() + 1);
		
		matchState(0, 0) = one;
		deletionState(0, 0) = zero;
		insertionState(0, 0) = zero;
		
		// Calculate first row
		if (seq2.size() > 0) {
			insertionState(0, 1) = insertionGap * insertionSpace;
			deletionState(0, 1) = matchState(0, 1) = zero;
			for (size_t j = 2; j <= seq2.size(); ++j) {
				insertionState(0, j) = insertionSpace * insertionState(0, j - 1);
				deletionState(0, j) = matchState(0, j) = zero;
			}
		}

		// Calculate all rows
		for (size_t i = 1; i <= seq1.size(); ++i) {
			for (size_t j = 0; j <= seq2.size(); ++j) {
				if (j == 0) {
					matchState(i, j) = insertionState(i, j) = zero;
				} else {
					matchState(i, j) = (matchState(i - 1, j - 1) +
										deletionState(i - 1, j - 1) +
										insertionState(i - 1, j - 1)) *
						matrix.getCharScore(seq1[i - 1], seq2[j - 1]);
					insertionState(i, j) =
						insertionSpace * (insertionGap *
										  (matchState(i, j - 1) +
										   deletionState(i, j - 1))
										  + insertionState(i, j - 1));
				}

				deletionState(i, j) =
					deletionSpace * (deletionGap * (matchState(i - 1, j) +
													insertionState(i - 1, j)) +
									 deletionState(i - 1, j));
			}
		}
	}
		
	
	template<typename SemiRing>
	void
	AffineGapNWPairwiseScorer<SemiRing>::
	scoreMatrixBackward(const std::string& seq1,
						const std::string& seq2,
						ScoreMatrix& matchState,
						ScoreMatrix& deletionState,
						ScoreMatrix& insertionState) const {
		const Score zero = semiRing.getZero();
		const Score one = semiRing.getMultiplicativeIdentity();

		size_t n = seq1.size();
		size_t m = seq2.size();
		
		// Initialize matrix
		matchState.resize(n + 1, m + 1);
		deletionState.resize(n + 1, m + 1);
		insertionState.resize(n + 1, m + 1);
		
		matchState(n, m) = deletionState(n, m) = insertionState(n, m) = one;

		// Calculate last row
		for (size_t j = 0; j < m; ++j) {
			size_t s = m - j; // source column
			size_t t = s - 1; // target column
			insertionState(n, t) = insertionState(n, s) * insertionSpace;
			matchState(n, t) = deletionState(n, t) = insertionState(n, t) * insertionGap;
		}

		// Calculate all rows
		for (size_t i = 1; i <= n; ++i) {
			size_t tRow = n - i; // target row
			size_t sRow = tRow + 1; // source row
			for (size_t j = 0; j <= m; ++j) {
				size_t tCol = m - j; // target column
				size_t sCol = tCol + 1; // source column

				deletionState(tRow, tCol) = deletionState(sRow, tCol) * deletionSpace;
				matchState(tRow, tCol) = insertionState(tRow, tCol)
					= deletionState(tRow, tCol) * deletionGap;
				
				if (tCol < m) {
					Score diag = matchState(sRow, sCol) * 
						matrix.getCharScore(seq1[tRow], seq2[tCol]);
					matchState(tRow, tCol) = insertionState(tRow, tCol) += diag;
					deletionState(tRow, tCol) += diag;

					Score horiz = insertionState(tRow, sCol) * insertionSpace;
					matchState(tRow, tCol) += horiz * insertionGap;
					deletionState(tRow, tCol) += horiz * insertionGap;
					insertionState(tRow, tCol) += horiz;
				}
			}
		}
	}

	template<typename SemiRing>
	void
	AffineGapNWPairwiseScorer<SemiRing>::
	scoreMatrix(const std::string& seq1,
				const std::string& seq2,
				ScoreMatrix& matchState,
				ScoreMatrix& deletionState,
				ScoreMatrix& insertionState) const {
		scoreMatrixForward(seq1, seq2, matchState, deletionState, insertionState);
		ScoreMatrix bMatchState, bDeletionState, bInsertionState;
		scoreMatrixBackward(seq1, seq2, bMatchState, bDeletionState, bInsertionState);
		matchState *= bMatchState;
		deletionState *= bDeletionState;
		insertionState *= bInsertionState;
	}
	
	template<typename SemiRing>
	typename SemiRing::Element
	AffineGapNWPairwiseScorer<SemiRing>::
	score(const std::string& seq1,
		  const std::string& seq2) const {
		ScoreList fromMatch, fromDeletion, fromInsertion;
		scoreLastRow(seq1, seq2, fromMatch, fromDeletion, fromInsertion, MATCH);
		return fromMatch.back() + fromDeletion.back() + fromInsertion.back();
	}
		
} }

#endif // __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISESCORER_HH__
