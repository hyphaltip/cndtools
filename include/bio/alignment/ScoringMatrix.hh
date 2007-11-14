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

#ifndef __BIO_ALIGNMENT_SCORINGMATRIX_HH__
#define __BIO_ALIGNMENT_SCORINGMATRIX_HH__

#include "bio/alphabet/Alphabet.hh"

namespace bio { namespace alignment {

	template<typename Score>
	class ScoringMatrix {
	public:
		virtual const alphabet::Alphabet& getAlphabet() const = 0;
		virtual Score getScore(char c1, char c2) const = 0;
		virtual Score getCharScore(char c1, char c2) const;
		virtual ~ScoringMatrix() {}
	};

	template<typename Score>
	Score ScoringMatrix<Score>::getCharScore(char c1, char c2) const {
		return getScore(getAlphabet().encode(c1), getAlphabet().encode(c2));
	}
	
} }

#endif // __BIO_ALIGNMENT_SCORINGMATRIX_HH__
