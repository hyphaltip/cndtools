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

#ifndef __AMINO_ACID_SCORING_MATRIX_HH__
#define __AMINO_ACID_SCORING_MATRIX_HH__

#include <istream>

#include "bio/alphabet/AminoAcid.hh"
#include "bio/alignment/AlphabetScoringMatrix.hh"

namespace bio { namespace alignment {

	template<typename Score>
	class AminoAcidScoringMatrix : public AlphabetScoringMatrix<Score> {
	public:
		AminoAcidScoringMatrix();
		AminoAcidScoringMatrix(Score match, Score mismatch);
		AminoAcidScoringMatrix(std::istream& strm);
	};

	template<typename Score>
	AminoAcidScoringMatrix<Score>::AminoAcidScoringMatrix()
		: AlphabetScoringMatrix<Score>(bio::alphabet::UnambiguousAminoAcid) {
	}

	template<typename Score>
	AminoAcidScoringMatrix<Score>::AminoAcidScoringMatrix(Score match, Score mismatch)
		: AlphabetScoringMatrix<Score>(bio::alphabet::UnambiguousAminoAcid, match, mismatch) {
	}

	template<typename Score>
	AminoAcidScoringMatrix<Score>::AminoAcidScoringMatrix(std::istream& strm)
		: AlphabetScoringMatrix<Score>(bio::alphabet::UnambiguousAminoAcid, strm) {
	}
	
} }

#endif // __AMINO_ACID_SCORING_MATRIX_HH__
