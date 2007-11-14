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

#ifndef __ALPHABET_SCORING_MATRIX_HH__
#define __ALPHABET_SCORING_MATRIX_HH__

#include <iosfwd>

#include "bio/alphabet/Alphabet.hh"
#include "bio/alignment/ScoringMatrix.hh"
#include "util/matrix.hh"

namespace bio { namespace alignment {

	template<typename Score>
	class AlphabetScoringMatrix : public ScoringMatrix<Score> {
	public:
		AlphabetScoringMatrix(const ScoringMatrix<Score>& matrix);
		AlphabetScoringMatrix(const alphabet::Alphabet& alphabet);
		AlphabetScoringMatrix(const alphabet::Alphabet& alphabet,
							  Score match,
							  Score mismatch);
		AlphabetScoringMatrix(const alphabet::Alphabet& alphabet,
							  std::istream& strm);

		AlphabetScoringMatrix& operator=(const ScoringMatrix<Score>& matrix);
		
		const alphabet::Alphabet& getAlphabet() const;
		
		Score getScore(char c1, char c2) const;
		Score getCharScore(char c1, char c2) const;

		void setScore(char c1, char c2, Score score);
		void setCharScore(char c1, char c2, Score score);
		void setScoreSymmetric(char c1, char c2, Score score);
		void setCharScoreSymmetric(char c1, char c2, Score score);
		void setMatchScore(Score score);
		void setMismatchScore(Score score);
		void readMatrix(std::istream& strm);

		template<typename S>
		friend std::ostream& operator<<(std::ostream& strm,
										const AlphabetScoringMatrix<S>& m);
		
	protected:
		const alphabet::Alphabet& alphabet;
		util::Matrix<Score> scores;
	};

	template<typename Score>
	AlphabetScoringMatrix<Score>::
	AlphabetScoringMatrix(const ScoringMatrix<Score>& matrix)
		: alphabet(matrix.getAlphabet()),
		  scores(alphabet.getSize(), alphabet.getSize()) {
		for (size_t i = 0; i < alphabet.getSize(); ++i) {
			for (size_t j = 0; j < alphabet.getSize(); ++j) {
				scores(i, j) = matrix.getScore(i, j);
			}
		}
	}

	template<typename Score>
	AlphabetScoringMatrix<Score>&
	AlphabetScoringMatrix<Score>::
	operator=(const ScoringMatrix<Score>& matrix) {
		for (size_t i = 0; i < alphabet.getSize(); ++i) {
			for (size_t j = 0; j < alphabet.getSize(); ++j) {
				scores(i, j) = matrix.getScore(i, j);
			}
		}
		return *this;
	}
	
	template<typename Score>
	AlphabetScoringMatrix<Score>::AlphabetScoringMatrix(const alphabet::Alphabet& alphabet)
		: alphabet(alphabet), scores(alphabet.getSize(), alphabet.getSize()) {
	}

	template<typename Score>
	AlphabetScoringMatrix<Score>::AlphabetScoringMatrix(const alphabet::Alphabet& alphabet,
														Score match,
														Score mismatch)
		: alphabet(alphabet), scores(alphabet.getSize(), alphabet.getSize()) {
		setMatchScore(match);
		setMismatchScore(mismatch);
	}

	template<typename Score>
	AlphabetScoringMatrix<Score>::AlphabetScoringMatrix(const alphabet::Alphabet& alphabet,
														std::istream& strm)
		: alphabet(alphabet), scores(alphabet.getSize(), alphabet.getSize()) {
		readMatrix(strm);
	}

	template<typename Score>
	const alphabet::Alphabet&	
	AlphabetScoringMatrix<Score>::getAlphabet() const { return alphabet; }
	
	template<typename Score>
	inline Score AlphabetScoringMatrix<Score>::getScore(char c1, char c2) const {
		return scores(c1, c2);
	}
	
	template<typename Score>
	inline Score AlphabetScoringMatrix<Score>::getCharScore(char c1, char c2) const {
		return scores(alphabet.encode(c1), alphabet.encode(c2));
	}
	
	template<typename Score>
	void AlphabetScoringMatrix<Score>::setScore(char c1, char c2, Score score) {
		scores(c1, c2) = score;
	}

	template<typename Score>
	void AlphabetScoringMatrix<Score>::setCharScore(char c1, char c2, Score score) {
		scores(alphabet.encode(c1), alphabet.encode(c2)) = score;
	}

	template<typename Score>
	void AlphabetScoringMatrix<Score>::setScoreSymmetric(char c1, char c2, Score score) {
		setScore(c1, c2, score);
		setScore(c2, c1, score);
	}
	
	template<typename Score>
	void AlphabetScoringMatrix<Score>::setCharScoreSymmetric(char c1, char c2, Score score) {
		setCharScore(c1, c2, score);
		setCharScore(c2, c1, score);
	}

	template<typename Score>
	void AlphabetScoringMatrix<Score>::setMatchScore(Score score) {
		for (size_t i = 0; i < alphabet.getSize(); ++i) {
			scores(i, i) = score;
		}
	}

	template<typename Score>
	void AlphabetScoringMatrix<Score>::setMismatchScore(Score score) {
		for (size_t i = 0; i < alphabet.getSize(); ++i) {
			for (size_t j = i + 1; j < alphabet.getSize(); ++j) {
				scores(i, j) = score;
				scores(j, i) = score;
			}
		}
	}

	template<typename Score>
	void AlphabetScoringMatrix<Score>::readMatrix(std::istream& strm) {
		for (size_t i = 0; i < alphabet.getSize(); ++i) {
			for (size_t j = 0; j < alphabet.getSize(); ++j) {
				strm >> scores(i, j);
			}
		}
	}

	template<typename Score>
	std::ostream& operator<<(std::ostream& strm,
							 const AlphabetScoringMatrix<Score>& m) {
		return strm << m.scores;
	}
	
} }

#endif // __ALPHABET_SCORING_MATRIX_HH__
