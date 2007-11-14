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

#ifndef __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISESUMMARIZER_HH__
#define __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISESUMMARIZER_HH__

#include "bio/alignment/AffineGapNWPairwiseScorer.hh"
#include "bio/alignment/AffineGapNWPairwiseAligner.hh"
#include "bio/alignment/AffineGapNWFastPairwiseAligner.hh"
#include "math/MaxPlus.hh"

namespace bio { namespace alignment {

	template<typename ScoreSemiRing>
	class SummarySemiRing {
	public:
		typedef typename ScoreSemiRing::Element Score;
		typedef std::vector<size_t> Summary;

		class Element {
		public:
			Score score;
			Summary summary;

			Element(Score score = 0, Summary summary = Summary());
			Element& operator=(const Element& x);
			Element& operator*=(const Element& x);
			Element& operator+=(const Element& x);
			Element operator*(const Element& x) const;
			Element operator+(const Element& x) const;
			bool operator<(const Element& x) const;
		};
		
		SummarySemiRing(const ScoreSemiRing& scoreSemiRing, size_t dims) :
			scoreSemiRing(scoreSemiRing), dims(dims) {
		}

		Element getZero() const {
			return Element(scoreSemiRing.getZero(),
						   Summary(dims));
		}

		Element getMultiplicativeIdentity() const {
			return Element(scoreSemiRing.getMultiplicativeIdentity(),
						   Summary(dims));
		}

		Element getUnitElement(Score score, size_t dim) {
			Summary summary(dims);
			summary[dim] = 1;
			return Element(score, summary);
		}
		
		ScoreSemiRing scoreSemiRing;
		size_t dims;
	};

	// I can't seem to get the code below to be recognized by the compiler
	//
	// 	template<typename ScoreSemiRing>
	// 	std::ostream& operator<<(std::ostream& strm, const SummarySemiRing<ScoreSemiRing>::Element& e) {
	// 		strm << '(' << e.score;
	// 		for (size_t i = 0; i < e.summary.size(); ++i) {
	// 			strm << ',' << ' ' << e.summary[i];
	// 		}
	// 		strm << ')';
	// 		return strm;
	// 	}

	template<typename ScoreSemiRing>
	SummarySemiRing<ScoreSemiRing>::Element::
	Element(Score score, Summary summary)
		: score(score), summary(summary) {
	}
	
	template<typename ScoreSemiRing>
	typename SummarySemiRing<ScoreSemiRing>::Element&
	SummarySemiRing<ScoreSemiRing>::Element::
	operator=(const Element& x) {
		if (summary.size() != x.summary.size()) {
			summary.resize(x.summary.size());
		}
		for (size_t i = 0; i < summary.size(); ++i) {
			summary[i] = x.summary[i];
		}
		score = x.score;
		return *this;
	}

	template<typename ScoreSemiRing>
	typename SummarySemiRing<ScoreSemiRing>::Element&
	SummarySemiRing<ScoreSemiRing>::Element::
	operator*=(const Element& x) {
		assert(summary.size() == x.summary.size());
		for (size_t i = 0; i < summary.size(); ++i) {
			summary[i] += x.summary[i];
		}
		score *= x.score;
		return *this;
	}

	template<typename ScoreSemiRing>
	typename SummarySemiRing<ScoreSemiRing>::Element&
	SummarySemiRing<ScoreSemiRing>::Element::
	operator+=(const Element& x) {
		if (*this < x) {
			*this = x;
		}
		return *this;
	}
	
	template<typename ScoreSemiRing>
	inline
	typename SummarySemiRing<ScoreSemiRing>::Element
	SummarySemiRing<ScoreSemiRing>::Element::
	operator*(const Element& x) const {
		return Element(*this) *= x;
	}
	
	template<typename ScoreSemiRing>
	inline
	typename SummarySemiRing<ScoreSemiRing>::Element
	SummarySemiRing<ScoreSemiRing>::Element::
	operator+(const Element& x) const {
		return Element(*this) += x;
	}
	
	template<typename ScoreSemiRing>
	bool
	SummarySemiRing<ScoreSemiRing>::Element::
	operator<(const Element& x) const {
		if (score < x.score) {
			return true;
		} else if (x.score < score) {
			return false;
		} else {
			assert(summary.size() == x.summary.size());
			for (size_t i = 0; i < summary.size(); ++i) {
				if (x.summary[i] < summary[i]) {
					return true;
				} else if (summary[i] < x.summary[i]) {
					return false;
				}
			}
			return false;
		}
	}	
	
	template<typename Score>
	class VarMappingScoringMatrix : public ScoringMatrix<Score> {
	public:
		VarMappingScoringMatrix(const ScoringMatrix<size_t>& matrix,
								const std::vector<Score>& scores)
			: matrix(matrix), scores(scores) {
		}
		
		Score getScore(char c1, char c2) const {
			return scores[matrix.getScore(c1, c2)];
		}

		const alphabet::Alphabet& getAlphabet() const { return matrix.getAlphabet(); }
		
	private:
		const ScoringMatrix<size_t>& matrix;
		std::vector<Score> scores;
	};

	template<typename NumType>
	class AffineGapNWPairwiseSummarizer {
	public:
		typedef typename math::MaxPlus<NumType> ScoringSemiRing;
		typedef SummarySemiRing<ScoringSemiRing> SemiRing;
		typedef typename SemiRing::Summary Summary;
		typedef typename SemiRing::Element Element;

		AffineGapNWPairwiseSummarizer(const alphabet::Alphabet& alphabet);
		
		AffineGapNWPairwiseSummarizer(const ScoringMatrix<size_t>& matrixVars,
									  size_t spaceVar,
									  size_t gapVar,
									  bool isLastVarZero = false);

		void setIndices(const ScoringMatrix<size_t>& matrixVars,
						size_t spaceVar,
						size_t gapVar,
						bool isLastVarZero = false);
		
		Summary summarize(const std::string& seq1,
						  const std::string& seq2,
						  bool lexMin = false,
						  bool linearMem = false) const;

		Summary summarize(const PairwiseAlignment& alignment) const;
		
		void setScores(const std::vector<NumType>& scores);
		
	protected:
		Summary summarizeAny(const std::string& seq1,
							 const std::string& seq2,
							 bool linearMem) const;
		
		Summary summarizeLexMin(const std::string& seq1,
								const std::string& seq2) const;

		AlphabetScoringMatrix<size_t> matrixVars;
		size_t spaceVar;
		size_t gapVar;
		bool isLastVarZero;
		std::vector<NumType> scores;
		mutable AffineGapNWFastPairwiseAligner<NumType> fastAligner;		
	};
	
	template<typename NumType>
	AffineGapNWPairwiseSummarizer<NumType>::
	AffineGapNWPairwiseSummarizer(const ScoringMatrix<size_t>& matrixVars,
								  size_t spaceVar,
								  size_t gapVar,
								  bool isLastVarZero)
		: matrixVars(matrixVars),
		  spaceVar(spaceVar),
		  gapVar(gapVar),
		  isLastVarZero(isLastVarZero),
		  scores(),
		  fastAligner(matrixVars.getAlphabet()) {
	}

	template<typename NumType>
	AffineGapNWPairwiseSummarizer<NumType>::
	AffineGapNWPairwiseSummarizer(const alphabet::Alphabet& alphabet)
		: matrixVars(alphabet),
		  fastAligner(alphabet) {
	}
	
	template<typename NumType>
	void
	AffineGapNWPairwiseSummarizer<NumType>::
	setIndices(const ScoringMatrix<size_t>& matrixVars,
			   size_t spaceVar,
			   size_t gapVar,
			   bool isLastVarZero) {
		this->matrixVars = matrixVars;
		this->spaceVar = spaceVar;
		this->gapVar = gapVar;
		this->isLastVarZero = isLastVarZero;
	}
	
	template<typename NumType>
	typename AffineGapNWPairwiseSummarizer<NumType>::Summary
	AffineGapNWPairwiseSummarizer<NumType>::
	summarize(const std::string& seq1,
			  const std::string& seq2,
			  bool lexMin,
			  bool linearMem) const {
		if (lexMin) {
			return summarizeLexMin(seq1, seq2);
		} else {
			return summarizeAny(seq1, seq2, linearMem);
		}
	}

	template<typename NumType>
	typename AffineGapNWPairwiseSummarizer<NumType>::Summary
	AffineGapNWPairwiseSummarizer<NumType>::
	summarizeAny(const std::string& seq1,
				 const std::string& seq2,
				 bool linearMem) const {
		std::vector<NumType> eltScores(scores);
		if (isLastVarZero) {
			eltScores.push_back(0);
		}
		VarMappingScoringMatrix<NumType> matrix(matrixVars, eltScores);
		if (linearMem) {
			AffineGapNWPairwiseAligner<NumType> aligner(matrix,
														eltScores[spaceVar],
														eltScores[gapVar]);
			return summarize(aligner.align(seq1, seq2));
		} else {
			fastAligner.setScores(matrix, eltScores[spaceVar], eltScores[gapVar]);
			return summarize(fastAligner.align(seq1, seq2));
		}
	}
		
	template<typename NumType>
	typename AffineGapNWPairwiseSummarizer<NumType>::Summary
	AffineGapNWPairwiseSummarizer<NumType>::
	summarizeLexMin(const std::string& seq1,
					const std::string& seq2) const {
		ScoringSemiRing scoringSemiRing;
		SemiRing summarySemiRing(scoringSemiRing, scores.size());
		std::vector<Element> eltScores;
		for (size_t i = 0; i < scores.size(); ++i) {
			eltScores.push_back(summarySemiRing.getUnitElement(scores[i], i));
		}
		if (isLastVarZero) {
			eltScores.push_back(summarySemiRing.getMultiplicativeIdentity());
		}
		VarMappingScoringMatrix<Element> matrix(matrixVars, eltScores);
		AffineGapNWPairwiseScorer<SemiRing>
			scorer(summarySemiRing,
				   matrix,
				   eltScores[spaceVar],
				   eltScores[gapVar]);
		return scorer.score(seq1, seq2).summary;
	}
	
	template<typename NumType>
	typename AffineGapNWPairwiseSummarizer<NumType>::Summary
	AffineGapNWPairwiseSummarizer<NumType>::
	summarize(const PairwiseAlignment& alignment) const {
		ScoringSemiRing scoringSemiRing;
		SemiRing summarySemiRing(scoringSemiRing, scores.size());
		std::vector<Element> eltScores;
		for (size_t i = 0; i < scores.size(); ++i) {
			eltScores.push_back(summarySemiRing.getUnitElement(scores[i], i));
		}
		if (isLastVarZero) {
			eltScores.push_back(summarySemiRing.getMultiplicativeIdentity());
		}
		VarMappingScoringMatrix<Element> matrix(matrixVars, eltScores);
		Element spaceScore = eltScores[spaceVar];
		Element gapScore = eltScores[gapVar];

		bool inGap1 = false;
		bool inGap2 = false;
		Element result = summarySemiRing.getMultiplicativeIdentity();
		for (size_t i = 0; i < alignment.length(); ++i) {
			if (alignment.seq1[i] == '-') {
				// Skip over positions that are gaps in both sequences
				// (useful for pairwise scoring of multiple
				// alignments)
				if (alignment.seq2[i] == '-') {
					continue;
				}
				result *= spaceScore;
				if (not inGap1) {
					result *= gapScore;
				}
				inGap1 = true;
				inGap2 = false;
			} else if (alignment.seq2[i] == '-') {
				result *= spaceScore;
				if (not inGap2) {
					result *= gapScore;
				}
				inGap1 = false;
				inGap2 = true;
			} else {
				result *= matrix.getCharScore(alignment.seq1[i], alignment.seq2[i]);
				inGap1 = false;
				inGap2 = false;
			}
		}

		return result.summary;
	}

	template<typename NumType>
	void
	AffineGapNWPairwiseSummarizer<NumType>::
	setScores(const std::vector<NumType>& scores) {
		this->scores = scores;
	}

} }

#endif // __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISESUMMARIZER_HH__
