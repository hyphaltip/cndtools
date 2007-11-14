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

#ifndef __AMBIGUOUS_DNA_SCORING_MATRIX_HH__
#define __AMBIGUOUS_DNA_SCORING_MATRIX_HH__

#include <istream>

#include "bio/alphabet/AmbiguousNucleotide.hh"
#include "bio/alignment/AlphabetScoringMatrix.hh"
#include "bio/alignment/DNAScoringMatrix.hh"

using bio::alphabet::AmbiguousDNA;

namespace bio { namespace alignment {

	template<typename Score>
	class AmbiguousDNAScoringMatrix : public AlphabetScoringMatrix<Score> {
	private:
		void init(const DNAScoringMatrix<Score>& dnaMatrix);
		
	public:
		AmbiguousDNAScoringMatrix();
		AmbiguousDNAScoringMatrix(const DNAScoringMatrix<Score>& dnaMatrix);
		AmbiguousDNAScoringMatrix(Score match, Score mismatch);
		
	};

	template<typename Score>
	void AmbiguousDNAScoringMatrix<Score>::init(const DNAScoringMatrix<Score>& dnaMatrix) {
		for (size_t i = 0; i < AmbiguousDNA.getSize(); ++i) {
			char ambiguousChar1 = AmbiguousDNA.decode(i);
			std::string chars1 = AmbiguousDNA.getPossibleChars(ambiguousChar1);
			for (size_t j = 0; j < AmbiguousDNA.getSize(); ++j) {
				char ambiguousChar2 = AmbiguousDNA.decode(j);
				std::string chars2 = AmbiguousDNA.getPossibleChars(ambiguousChar2);
				Score total = 0;
				typedef std::string::const_iterator CharIt;
				for (CharIt c1 = chars1.begin(); c1 != chars1.end(); ++c1) {
					for (CharIt c2 = chars2.begin(); c2 != chars2.end(); ++c2) {
						total += dnaMatrix.getCharScore(*c1, *c2);
					}
				}
				
				setCharScore(ambiguousChar1, ambiguousChar2,
							 total / static_cast<int>(chars1.size() * chars2.size()));
			}
		}
	}
	
	template<typename Score>
	AmbiguousDNAScoringMatrix<Score>::AmbiguousDNAScoringMatrix()
		: AlphabetScoringMatrix<Score>(AmbiguousDNA) {
	}

	template<typename Score>
	AmbiguousDNAScoringMatrix<Score>::AmbiguousDNAScoringMatrix(const DNAScoringMatrix<Score>& dnaMatrix)
		: AlphabetScoringMatrix<Score>(AmbiguousDNA) {
		init(dnaMatrix);
	}
	
	template<typename Score>
	AmbiguousDNAScoringMatrix<Score>::AmbiguousDNAScoringMatrix(Score match, Score mismatch)
		: AlphabetScoringMatrix<Score>(AmbiguousDNA) {
		DNAScoringMatrix<Score> dnaMatrix(match, mismatch);
		init(dnaMatrix);
	}

	
} }

#endif // __AMBIGUOUS_DNA_SCORING_MATRIX_HH__
