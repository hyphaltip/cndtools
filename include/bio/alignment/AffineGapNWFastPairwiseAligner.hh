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

#ifndef __BIO_ALIGNMENT_AFFINEGAPNWFASTPAIRWISEALIGNER_HH__
#define __BIO_ALIGNMENT_AFFINEGAPNWFASTPAIRWISEALIGNER_HH__

#include "bio/alignment/PairwiseAligner.hh"
#include "bio/alignment/AlphabetScoringMatrix.hh"

namespace bio { namespace alignment {

	template<typename NumType>	
    class AffineGapNWFastPairwiseAligner : public PairwiseAligner {
	public:
		AffineGapNWFastPairwiseAligner(const alphabet::Alphabet& alphabet);
		
		AffineGapNWFastPairwiseAligner(const ScoringMatrix<NumType>& matrix,
									   const NumType& space,
									   const NumType& gap);
		
		PairwiseAlignment align(const std::string& seq1,
								const std::string& seq2) const;

		void setScores(const ScoringMatrix<NumType>& matrix,
					   const NumType& space,
					   const NumType& gap);
		
	private:
		static PairwiseAlignment
		statesToAlignment(const std::vector<int>& indices,
						  const std::string& seq1,
						  const std::string& seq2);

		static int argmax(NumType h, NumType d, NumType i);

		static const int H;
		static const int D;
		static const int I;

		AlphabetScoringMatrix<NumType> matrix;
		NumType space;
		NumType gap;
		mutable util::Matrix<NumType> hState;
		mutable util::Matrix<NumType> dState;
		mutable util::Matrix<NumType> iState;
	};

	template<typename NumType>	
    const int AffineGapNWFastPairwiseAligner<NumType>::H = 0;
	template<typename NumType>	
    const int AffineGapNWFastPairwiseAligner<NumType>::D = 1;
	template<typename NumType>	
    const int AffineGapNWFastPairwiseAligner<NumType>::I = 2;

	template<typename NumType>
	AffineGapNWFastPairwiseAligner<NumType>::
	AffineGapNWFastPairwiseAligner(const alphabet::Alphabet& alphabet) 
		: matrix(alphabet),
		  space(),
		  gap() {
	}
	
	template<typename NumType>
	AffineGapNWFastPairwiseAligner<NumType>::
	AffineGapNWFastPairwiseAligner(const ScoringMatrix<NumType>& matrix,
								   const NumType& space,
								   const NumType& gap)
		: matrix(matrix),
		  space(space),
		  gap(gap) {
	}

	template<typename NumType>
	void
	AffineGapNWFastPairwiseAligner<NumType>::
	setScores(const ScoringMatrix<NumType>& matrix,
			  const NumType& space,
			  const NumType& gap) {
		this->matrix = matrix;
		this->space = space;
		this->gap = gap;
	}
	
	template<typename NumType>
	int
	AffineGapNWFastPairwiseAligner<NumType>::
	argmax(NumType h, NumType d, NumType i) {
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
	AffineGapNWFastPairwiseAligner<NumType>::
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
	AffineGapNWFastPairwiseAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2) const {
		using std::max;

		std::string eSeq1(matrix.getAlphabet().encode(seq1));
		std::string eSeq2(matrix.getAlphabet().encode(seq2));
		
		size_t n = seq1.size();
		size_t m = seq2.size();

		hState.resize(n + 1, m + 1);
		dState.resize(n + 1, m + 1);
		iState.resize(n + 1, m + 1);
		
		hState(n, m) = dState(n, m) = iState(n, m) = 0;

		// Calculate last row
		for (size_t j = 0; j < m; ++j) {
			size_t s = m - j; // source column
			size_t t = s - 1; // target column
			iState(n, t) = iState(n, s) + space;
			hState(n, t) = dState(n, t) = iState(n, t) + gap;
		}

		// Calculate all rows
		for (size_t i = 1; i <= n; ++i) {
			size_t tRow = n - i; // target row
			size_t sRow = tRow + 1; // source row
			for (size_t j = 0; j <= m; ++j) {
				size_t tCol = m - j; // target column
				size_t sCol = tCol + 1; // source column

				NumType dVal = dState(sRow, tCol) + space;
				NumType hVal = dVal + gap;
				NumType iVal = hVal;

				if (tCol < m) {
					NumType diag = hState(sRow, sCol)
						+ matrix.getScore(eSeq1[tRow], eSeq2[tCol]);
					hVal = iVal = max(diag, hVal);
					dVal = max(diag, dVal);
					
					NumType horiz = iState(tRow, sCol) + space;
					hVal = max(horiz + gap, hVal);
					dVal = max(horiz + gap, dVal);
					iVal = max(horiz, iVal);
				}

				hState(tRow, tCol) = hVal;
				dState(tRow, tCol) = dVal;
				iState(tRow, tCol) = iVal;				
			}
		}
	
		std::vector<int> traceback;

		size_t i = 0, j = 0;
		
		int currState = H;

		while (i < n or j < m) {
			if (i == n) {
				currState = I;
			} else if (j == m) {
				currState = D;
			} else {
				NumType hScore = hState(i + 1, j + 1)
					+ matrix.getScore(eSeq1[i], eSeq2[j]);
				NumType dScore = dState(i + 1, j) + space;
				NumType iScore = iState(i, j + 1) + space;

				if (currState != D) { dScore += gap; }
				if (currState != I) { iScore += gap; }

				currState = argmax(hScore, dScore, iScore);
			}

			if (currState == H) {
				++i;
				++j;
			} else if (currState == D) {
				++i;
			} else if (currState == I) {
				++j;
			} else {
				throw std::runtime_error("AffineGapNWFastPairwiseAligner"
										 "::align: Bad state");
			}

			traceback.push_back(currState);
		}
		
		return statesToAlignment(traceback, seq1, seq2);
	}
	
} }

#endif // __BIO_ALIGNMENT_AFFINEGAPNWFASTPAIRWISEALIGNER_HH__
