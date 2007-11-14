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

#ifndef __DNA_SCORING_MATRIX_HH__
#define __DNA_SCORING_MATRIX_HH__

#include <iosfwd>

#include "bio/alphabet/Nucleotide.hh"
#include "bio/alignment/AlphabetScoringMatrix.hh"

namespace bio { namespace alignment {

	template<typename Score>
	class DNAScoringMatrix : public AlphabetScoringMatrix<Score> {
	public:
		DNAScoringMatrix();
		DNAScoringMatrix(Score match, Score mismatch);
		DNAScoringMatrix(Score match,
						 Score transition,
						 Score transversion);
		DNAScoringMatrix(std::istream& strm);

		void setTransversionScore(Score s);
		void setTransitionScore(Score s);
	};

	template<typename Score>
	DNAScoringMatrix<Score>::DNAScoringMatrix()
		: AlphabetScoringMatrix<Score>(bio::alphabet::DNA) {
	}

	template<typename Score>
	DNAScoringMatrix<Score>::DNAScoringMatrix(Score match, Score mismatch)
		: AlphabetScoringMatrix<Score>(bio::alphabet::DNA, match, mismatch) {
	}

	template<typename Score>
	DNAScoringMatrix<Score>::DNAScoringMatrix(std::istream& strm)
		: AlphabetScoringMatrix<Score>(bio::alphabet::DNA, strm) {
	}

	template<typename Score>
	DNAScoringMatrix<Score>::DNAScoringMatrix(Score match,
											  Score transition,
											  Score transversion)
		: AlphabetScoringMatrix<Score>(bio::alphabet::DNA) {
		setMatchScore(match);
		setTransitionScore(transition);
		setTransversionScore(transversion);
	}

	template<typename Score>
	void
	DNAScoringMatrix<Score>::setTransversionScore(Score s) {
		using alphabet::DNA_PURINES;
		using alphabet::DNA_PYRIMIDINES;
		for (size_t i = 0; i < DNA_PURINES.getSize(); ++i) {
			for (size_t j = 0; j < DNA_PYRIMIDINES.getSize(); ++j) {
				setCharScoreSymmetric(DNA_PURINES.decode(i),
									  DNA_PYRIMIDINES.decode(j),
									  s);
			}
		}
	}

	template<typename Score>
	void
	DNAScoringMatrix<Score>::setTransitionScore(Score s) {
		using alphabet::DNA_PURINES;
		using alphabet::DNA_PYRIMIDINES;
		for (size_t i = 0; i < DNA_PURINES.getSize(); ++i) {
			for (size_t j = i + 1; j < DNA_PURINES.getSize(); ++j) {
				setCharScoreSymmetric(DNA_PURINES.decode(i),
									  DNA_PURINES.decode(j),
									  s);
			}
		}

		for (size_t i = 0; i < DNA_PYRIMIDINES.getSize(); ++i) {
			for (size_t j = i + 1; j < DNA_PYRIMIDINES.getSize(); ++j) {
				setCharScoreSymmetric(DNA_PYRIMIDINES.decode(i),
									  DNA_PYRIMIDINES.decode(j),
									  s);
			}
		}
	}
		
} }

#endif // __DNA_SCORING_MATRIX_HH__
