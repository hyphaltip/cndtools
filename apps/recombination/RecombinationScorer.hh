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

#include <vector>

#include "util/matrix.hh"
#include "bio/alignment/AlphabetScoringMatrix.hh"
#include "bio/alignment/MultipleAlignment.hh"
using bio::alignment::AlphabetScoringMatrix;
using bio::alignment::MultipleAlignment;

template<typename SemiRing>
class RecombinationScorer {
public:
	typedef typename SemiRing::Element Score;
	
	RecombinationScorer(const SemiRing& semiRing,
						const std::vector<Score>& initialScores,
						const AlphabetScoringMatrix<Score>& emissionScores,
						const util::Matrix<Score>& transitionScores);

	Score score(const MultipleAlignment& ma,
				size_t recombinant = 0) const;

private:
	SemiRing semiRing;
	std::vector<Score> initialScores;
	AlphabetScoringMatrix<Score> emissionScores;
	util::Matrix<Score> transitionScores;
};

template<typename SemiRing>
RecombinationScorer<SemiRing>::
RecombinationScorer(const SemiRing& semiRing,
					const std::vector<Score>& initialScores,
					const AlphabetScoringMatrix<Score>& emissionScores,
					const util::Matrix<Score>& transitionScores)
	: semiRing(semiRing),
	  initialScores(initialScores),
	  emissionScores(emissionScores),
	  transitionScores(transitionScores) {
}
	
template<typename SemiRing>
typename SemiRing::Element
RecombinationScorer<SemiRing>::
score(const MultipleAlignment& ma,
	  size_t recombinant) const {
	const Score zero = semiRing.getZero();
	const Score one = semiRing.getMultiplicativeIdentity();

	if (ma.getNumCols() == 0) {
		return zero;
	}

	std::vector<Score> prev(ma.getNumSeqs(), zero);
	std::vector<Score> next(ma.getNumSeqs(), zero);

	// Calculate first column
	for (size_t i = 0; i < ma.getNumSeqs(); ++i) {
		if (i == recombinant) { continue; }
		next[i] = emissionScores.getCharScore(ma.getChar(recombinant, 0),
											  ma.getChar(i, 0));
		next[i] *= initialScores[i];
	}
	
	// Calculate the rest of the columns
	for (size_t t = 1; t < ma.getNumCols(); ++t) {
		
		next.swap(prev);
		
		for (size_t i = 0; i < ma.getNumSeqs(); ++i) {
			if (i == recombinant) { continue; }

			next[i] = zero;
			
			for (size_t j = 0; j < ma.getNumSeqs(); ++j) {
				if (j == recombinant) { continue; }
				next[i] += (prev[j] * transitionScores(j, i));
			}
			
			next[i] *= emissionScores.getCharScore(ma.getChar(recombinant, t),
												   ma.getChar(i, t));
		}
	}
	
	Score result = zero;
	for (size_t i = 0; i < ma.getNumSeqs(); ++i) {
		if (i == recombinant) { continue; }
		result += next[i];
	}
	
	return result;
}
