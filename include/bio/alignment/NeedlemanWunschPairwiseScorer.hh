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

#ifndef __NEEDLEMAN_WUNSCH_PAIRWISE_SCORER_HH__
#define __NEEDLEMAN_WUNSCH_PAIRWISE_SCORER_HH__

#include <vector>

#include "bio/alignment/PairwiseSequenceScorer.hh"
#include "bio/alignment/ScoringMatrix.hh"
#include "util/stl.hh"
#include "util/matrix.hh"

namespace bio { namespace alignment {

	using util::stl::print_elements;
	using util::Matrix;
	
	template<typename SemiRing>
	class NeedlemanWunschPairwiseScorer
		: PairwiseSequenceScorer<typename SemiRing::Element> {
	public:
		typedef typename SemiRing::Element Score;
		
		NeedlemanWunschPairwiseScorer(const ScoringMatrix<Score>& matrix,
									  const Score& space);
		
		Score score(const std::string& seq1,
					const std::string& seq2) const;

	protected:
		typedef typename std::vector<Score> ScoreList;
		typedef typename util::Matrix<Score> ScoreMatrix;
		
		void scoreLastRow(const std::string& seq1,
						  const std::string& seq2,
						  ScoreList& row) const;

		void scoreMatrix(const std::string& seq1,
						 const std::string& seq2,
						 ScoreMatrix& m) const;
		
		const ScoringMatrix<Score>& matrix;
		Score space;
	};

	template<typename SemiRing>
	NeedlemanWunschPairwiseScorer<SemiRing>::
	NeedlemanWunschPairwiseScorer(const ScoringMatrix<typename SemiRing::Element>& matrix,
								  const typename SemiRing::Element& space)
		: matrix(matrix), space(space) {
	}
		
	template<typename SemiRing>
	void
	NeedlemanWunschPairwiseScorer<SemiRing>::
	scoreLastRow(const std::string& seq1,
				 const std::string& seq2,
				 NeedlemanWunschPairwiseScorer<SemiRing>::ScoreList& row) const {
		const Score one = SemiRing::multiplicativeIdentity;
		
		// Initialize row
		row.resize(seq2.size() + 1);
		row[0] = one;

		// Calculate first row
		for (size_t j = 1; j <= seq2.size(); ++j) {
			row[j] = space * row[j - 1];
		}

		// 		std::cerr << "First row:" << std::endl;
		// 		print_elements(std::cerr, row.begin(), row.end());
		// 		std::cerr << '\n';
		
		// Calculate remaining rows
		Score up, diag;
		for (size_t i = 1; i <= seq1.size(); ++i) {
			diag = row[0];
			row[0] *= space;
			for (size_t j = 1; j <= seq2.size(); ++j) {
				up = row[j];

				row[j] = matrix.getCharScore(seq1[i - 1], seq2[j - 1]) * diag
					+ space * (up + row[j - 1]);

				diag = up;
			}
			
			// 			std::cerr << "Row: " << i << std::endl;
			// 			print_elements(std::cerr, row.begin(), row.end());
			// 			std::cerr << '\n';
		}
	}

	template<typename SemiRing>
	void
	NeedlemanWunschPairwiseScorer<SemiRing>::
	scoreMatrix(const std::string& seq1,
				const std::string& seq2,
				Matrix<typename SemiRing::Element>& m) const {
		const Score one = SemiRing::multiplicativeIdentity;
		
		// Initialize matrix
		m.resize(seq1.size() + 1, seq2.size() + 1);
	    m(0, 0) = one;

		// Calculate first row
		for (size_t j = 1; j <= seq2.size(); ++j) {
			m(0, j) = space * m(0, j - 1);
		}
		
		// Calculate remaining rows
		for (size_t i = 1; i <= seq1.size(); ++i) {
			m(i, 0) = space * m(i - 1, 0);
			for (size_t j = 1; j <= seq2.size(); ++j) {
				m(i, j) =
					matrix.getCharScore(seq1[i - 1], seq2[j - 1]) * m(i - 1, j - 1)
					+ space * (m(i - 1, j) + m(i, j - 1));
			}
		}
	}
	
	template<typename SemiRing>
	typename SemiRing::Element
	NeedlemanWunschPairwiseScorer<SemiRing>::
	score(const std::string& seq1,
		  const std::string& seq2) const {
		ScoreList row;
		scoreLastRow(seq1, seq2, row);
		return row.back();
	}
	
} }

#endif // __NEEDLEMAN_WUNSCH_PAIRWISE_SCORER_HH__
