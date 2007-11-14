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

#ifndef __BIO_ALIGNMENT_AFFINEGAPNWFASTPAIRWISEALIGNER2_HH__
#define __BIO_ALIGNMENT_AFFINEGAPNWFASTPAIRWISEALIGNER2_HH__

#include "bio/alignment/AffineGapNWPairwiseScorer.hh"
#include "bio/alignment/PairwiseAlignment.hh"
#include "boost/shared_ptr.hpp"
using boost::shared_ptr;

namespace bio { namespace alignment {

	struct Path;
	typedef shared_ptr<Path> PathPtr;

	struct Path {
		size_t index;
		PathPtr last;

		Path(size_t index = 0, PathPtr last = PathPtr())
			: index(index), last(last) {}

		Path(PathPtr head, PathPtr tail) {
			index = head->index;
			if (head->last.get() == NULL) {
				last = tail;
			} else {
				last = PathPtr(new Path(head->last, tail));
			}
		}
		
	};

	void get_path_indices(PathPtr path, std::vector<size_t>& indices) {
		while (path != NULL) {
			indices.push_back(path->index);
			path = path->last;
		}
		std::reverse(indices.begin(), indices.end());
	}

	template<typename Score>
	struct ScoredPath {
		Score score;
		PathPtr path;

		ScoredPath(Score score = Score(), PathPtr path = PathPtr())
			: score(score), path(path) {}

		ScoredPath& operator+=(const ScoredPath& x) {
			if (score < x.score) { *this = x; }
			return *this;
		}

		ScoredPath& operator*=(const ScoredPath& x) {
			score *= x.score;
			if (x.path != NULL) {
				if (path != NULL) {
					path = PathPtr(new Path(x.path, path));
				} else {
					path = x.path;
				}
			}
			return *this;
		}

		ScoredPath operator*(const ScoredPath& x) {
			return ScoredPath(*this) *= x;
		}

		ScoredPath operator+(const ScoredPath& x) {
			return ScoredPath(*this) += x;
		}
	};

	template<typename ScoringSemiRing>
	class BestPathSemiRing {
	public:
		typedef ScoredPath<typename ScoringSemiRing::Element> Element;

		BestPathSemiRing(const ScoringSemiRing& scoringSemiRing)
			: scoringSemiRing(scoringSemiRing) {
		}
	
		Element getZero() const {
			return Element(scoringSemiRing.getZero());
		}
	
		Element getMultiplicativeIdentity() const {
			return Element(scoringSemiRing.getMultiplicativeIdentity());
		}

	private:
		ScoringSemiRing scoringSemiRing;
	};

	template<typename NumType, typename Score>
    class PathScoringMatrix : public ScoringMatrix<Score> {
	public:
		PathScoringMatrix(const ScoringMatrix<NumType>& matrix,
						  const PathPtr path)
			: matrix(matrix), path(path) {
		}
		
		Score getScore(char c1, char c2) const {
			return Score(matrix.getScore(c1, c2), path);
		}

		const alphabet::Alphabet& getAlphabet() const {
			return matrix.getAlphabet();
		}

	private:
		const ScoringMatrix<NumType>& matrix;
		const PathPtr path;
    };
	
	template<typename NumType>	
    class AffineGapNWFastPairwiseAligner2 : public PairwiseAligner {
	public:
		AffineGapNWFastPairwiseAligner2(const ScoringMatrix<NumType>& matrix,
										const NumType& space,
										const NumType& gap);
		
		PairwiseAlignment align(const std::string& seq1,
								const std::string& seq2) const;
		
	private:
		typedef math::MaxPlus<NumType> NumSemiRing;
		typedef BestPathSemiRing<NumSemiRing> PathSemiRing;
		typedef typename PathSemiRing::Element Score;
		typedef AffineGapNWPairwiseScorer<PathSemiRing> Scorer;

		static const size_t H;
		static const size_t I;
		static const size_t D;

		static PairwiseAlignment
		indicesToAlignment(const std::vector<size_t>& indices,
						   const std::string& seq1,
						   const std::string& seq2);
		
		PathScoringMatrix<NumType, Score> matrix;
		Score insertionSpace;
		Score deletionSpace;
		Score insertionGap;
		Score deletionGap;
	};

	template<typename NumType>	
    const size_t AffineGapNWFastPairwiseAligner2<NumType>::H = 0;
	template<typename NumType>	
    const size_t AffineGapNWFastPairwiseAligner2<NumType>::I = 1;
	template<typename NumType>	
    const size_t AffineGapNWFastPairwiseAligner2<NumType>::D = 2;

	template<typename NumType>
	AffineGapNWFastPairwiseAligner2<NumType>::
	AffineGapNWFastPairwiseAligner2(const ScoringMatrix<NumType>& matrix,
								   const NumType& space,
								   const NumType& gap)
		: matrix(matrix, PathPtr(new Path(H))),
		  insertionSpace(space, PathPtr(new Path(I))),
		  deletionSpace(space, PathPtr(new Path(D))),
		  insertionGap(gap),
		  deletionGap(gap) {
	}
	
	template<typename NumType>
	PairwiseAlignment
	AffineGapNWFastPairwiseAligner2<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2) const {
		NumSemiRing numSemiRing;
		PathSemiRing pathSemiRing(numSemiRing);
		Scorer scorer(pathSemiRing,
					  matrix,
					  insertionSpace,
					  deletionSpace,
					  insertionGap,
					  deletionGap);
		Score s = scorer.score(seq1, seq2);
		std::vector<size_t> indices;
		get_path_indices(s.path, indices);
		return indicesToAlignment(indices, seq1, seq2);
	}

	template<typename NumType>
	PairwiseAlignment
	AffineGapNWFastPairwiseAligner2<NumType>::
	indicesToAlignment(const std::vector<size_t>& indices,
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
	
	
	
} }

#endif // __BIO_ALIGNMENT_AFFINEGAPNWFASTPAIRWISEALIGNER2_HH__
