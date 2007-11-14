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

#ifndef __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISEALIGNER_HH__
#define __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISEALIGNER_HH__

#include "bio/alignment/PairwiseAligner.hh"
#include "bio/alignment/AffineGapNWPairwiseScorer.hh"
#include "bio/alignment/ConvertingScoringMatrix.hh"
#include "math/MaxPlus.hh"
#include "util/stl.hh"

namespace bio { namespace alignment {

	template<typename NumType>	
    class AffineGapNWPairwiseAligner : public PairwiseAligner {
	public:
		AffineGapNWPairwiseAligner(const ScoringMatrix<NumType>& matrix,
								   const NumType& space,
								   const NumType& gap);
		
		PairwiseAlignment align(const std::string& seq1,
								const std::string& seq2) const;

		void align(const std::string& seq1,
				   const std::string& seq2,
				   std::vector<PairwiseAlignment>& paths) const;
				
	private:
		typedef math::MaxPlus<NumType> NumSemiRing;
		typedef AffineGapNWPairwiseScorer<NumSemiRing> Scorer;
  		typedef typename Scorer::Score Score;
  		typedef typename Scorer::ScoreList ScoreList;
		
		PairwiseAlignment align(const std::string& seq1,
								const std::string& seq2,
								size_t startState,
								size_t endState) const;

		PairwiseAlignment align(char c1,
								char c2,
								size_t startState,
								size_t endState) const;

		void align(const std::string& seq1,
				   const std::string& seq2,
				   std::vector<PairwiseAlignment>& paths,
				   size_t startState,
				   size_t endState) const;

		void align(char c1,
				   char c2,
				   std::vector<PairwiseAlignment>& paths,
				   size_t startState,
				   size_t endState) const;

		void align(const std::string& seq1,
				   const std::string& seq2,
				   std::vector<PairwiseAlignment>& paths,
				   size_t startState,
				   size_t endState,
				   size_t midState,
				   std::string::size_type seq1Prefix,
				   std::string::size_type seq2Prefix) const;
		
		void combineScores(ScoreList& l1, ScoreList& l2) const;
		size_t flipState(const size_t state) const;

		ConvertingScoringMatrix<NumType, Score> convertingMatrix;
		Score space;
		Score gap;
		Scorer scorer;
	};

	template<typename NumType>
	AffineGapNWPairwiseAligner<NumType>::
	AffineGapNWPairwiseAligner(const ScoringMatrix<NumType>& matrix,
							   const NumType& space,
							   const NumType& gap)
		: convertingMatrix(matrix),
		  space(space), gap(gap),
		  scorer(NumSemiRing(), convertingMatrix, space, gap) {
	}
	
	template<typename NumType>
	void
	AffineGapNWPairwiseAligner<NumType>::
	combineScores(ScoreList& l1, ScoreList& l2) const {
		std::transform(l1.begin(), l1.end(), l2.begin(), l1.begin(),
					   std::multiplies<Score>());
	}

	template<typename NumType>
	PairwiseAlignment
	AffineGapNWPairwiseAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2) const {
		return align(seq1, seq2, Scorer::MATCH, Scorer::MATCH | Scorer::GAP1 | Scorer::GAP2);
	}

	template<typename NumType>
	void
	AffineGapNWPairwiseAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2,
		  std::vector<PairwiseAlignment>& paths) const {
		align(seq1, seq2, paths, Scorer::MATCH, Scorer::MATCH | Scorer::GAP1 | Scorer::GAP2);
	}
	
	template<typename NumType>
	size_t
	AffineGapNWPairwiseAligner<NumType>::
	flipState(const size_t state) const {
		return ((state & Scorer::MATCH ? Scorer::MATCH : 0) |
				(state & Scorer::GAP1 ? Scorer::GAP2 : 0) |
				(state & Scorer::GAP2 ? Scorer::GAP1 : 0));
	}
	
	template<typename NumType>
	PairwiseAlignment
	AffineGapNWPairwiseAligner<NumType>::
	align(char c1,
		  char c2,
		  size_t startState,
		  size_t endState) const {
		NumSemiRing semiRing;
		const Score zero = semiRing.getZero();
		
		Score
			gap1FirstPathValue = zero,
			gap2FirstPathValue = zero,
			matchPathValue = zero;

		if (endState & Scorer::MATCH) {
			matchPathValue = convertingMatrix.getCharScore(c1, c2);
		}
		if (endState & Scorer::GAP1) {
			gap2FirstPathValue = space * gap * space;
			if (not (startState & Scorer::GAP2)) {
				gap2FirstPathValue *= gap;
			}
		}
		if (endState & Scorer::GAP2) {
			gap1FirstPathValue = space * gap * space;
			if (not (startState & Scorer::GAP1)) {
				gap1FirstPathValue *= gap;
			}
		}
		if (gap1FirstPathValue <= matchPathValue and
			gap2FirstPathValue <= matchPathValue) {
			return PairwiseAlignment(std::string(1, c1),
									 std::string(1, c2));
		} else if (gap2FirstPathValue <= gap1FirstPathValue) {
			return PairwiseAlignment(c1 + std::string("-"),
									 std::string("-") + c2);
		} else {
			return PairwiseAlignment(std::string("-") + c1,
									 c2 + std::string("-"));
		}
	}

	template<typename NumType>
	void
	AffineGapNWPairwiseAligner<NumType>::
	align(char c1,
		  char c2,
		  std::vector<PairwiseAlignment>& paths,
		  size_t startState,
		  size_t endState) const {
		NumSemiRing semiRing;
		const Score zero = semiRing.getZero();
		
		Score
			gap1FirstPathValue = zero,
			gap2FirstPathValue = zero,
			matchPathValue = zero;

		if (endState & Scorer::MATCH) {
			matchPathValue = convertingMatrix.getCharScore(c1, c2);
		}
		if (endState & Scorer::GAP1) {
			gap2FirstPathValue = space * gap * space;
			if (not (startState & Scorer::GAP2)) {
				gap2FirstPathValue *= gap;
			}
		}
		if (endState & Scorer::GAP2) {
			gap1FirstPathValue = space * gap * space;
			if (not (startState & Scorer::GAP1)) {
				gap1FirstPathValue *= gap;
			}
		}

		Score maxScore =
			matchPathValue + gap1FirstPathValue + gap2FirstPathValue;
		
		if (matchPathValue == maxScore) {
			paths.push_back(PairwiseAlignment(std::string(1, c1),
											  std::string(1, c2)));
		}
		if (gap1FirstPathValue == maxScore) {
			paths.push_back(PairwiseAlignment(c1 + std::string("-"),
											  std::string("-") + c2));
		}
		if (gap2FirstPathValue == maxScore) {
			paths.push_back(PairwiseAlignment(std::string("-") + c1,
											  c2 + std::string("-")));
		}
	}
	
	template<typename NumType>
	PairwiseAlignment
	AffineGapNWPairwiseAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2,
		  size_t startState,
		  size_t endState) const {
		if (seq1.size() == 0) {
			return PairwiseAlignment(std::string(seq2.size(), '-'), seq2);
		} else if (seq2.size() == 0) {
			return PairwiseAlignment(seq1, std::string(seq1.size(), '-'));
		} else if (seq1.size() == 1 and seq2.size() == 1) {
			return align(seq1[0], seq2[0], startState, endState);
		} else if (seq1.size() == 1) {
			PairwiseAlignment flipped =
				align(seq2, seq1, flipState(startState), flipState(endState));
			flipped.flip();
			return flipped;
		} else {
			// Calculate middle of first sequence
			std::string::size_type middle = seq1.size() / 2;
			std::string seq1Top = seq1.substr(0, middle);
			std::string seq1Bot = seq1.substr(middle);

			// Calculate forward values to middle of first sequence
			ScoreList forwardMatch, forwardGap1, forwardGap2;
			scorer.scoreLastRow(seq1Top, seq2,
								forwardMatch, forwardGap1, forwardGap2, startState);

			// Calculate backward values to middle of first sequence
			ScoreList backwardMatch, backwardGap1, backwardGap2;
			scorer.scoreFirstRow(seq1Bot, seq2,
								 backwardMatch, backwardGap1, backwardGap2, endState);
			
			// Combine forward and backward scores
			combineScores(forwardMatch, backwardMatch);
			combineScores(forwardGap1, backwardGap1);
			combineScores(forwardGap2, backwardGap2);

			// Find maximum positions
			typename ScoreList::iterator maxMatch, maxGap1;
			maxMatch = std::max_element(forwardMatch.begin(), forwardMatch.end());
			maxGap1 = std::max_element(forwardGap1.begin(), forwardGap1.end());

			Score maxScore = *maxMatch + *maxGap1;

			size_t middleState;
			size_t seq2Prefix;
			if (*maxMatch == maxScore) {
				middleState = Scorer::MATCH;
				seq2Prefix = maxMatch - forwardMatch.begin();
			} else {
				middleState = Scorer::GAP1;
				seq2Prefix = maxGap1 - forwardGap1.begin();
			}

			std::string seq2Left = seq2.substr(0, seq2Prefix);
			std::string seq2Right = seq2.substr(seq2Prefix);
			
			// Combine best alignment of top with best alignment of bottom
			PairwiseAlignment pa =
				align(seq1Top, seq2Left, startState, middleState) +
				align(seq1Bot, seq2Right, middleState, endState);
			return pa;
		}
	}

	template<typename NumType>
	void
	AffineGapNWPairwiseAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2,
		  std::vector<PairwiseAlignment>& paths,
		  size_t startState,
		  size_t endState,
		  size_t midState,
		  std::string::size_type seq1Prefix,
		  std::string::size_type seq2Prefix) const {

		std::string seq1Top = seq1.substr(0, seq1Prefix);
		std::string seq1Bot = seq1.substr(seq1Prefix);
		std::string seq2Left = seq2.substr(0, seq2Prefix);
		std::string seq2Right = seq2.substr(seq2Prefix);
		
		// Calculate top subalignments
		AlignmentList topAligns;
		align(seq1Top, seq2Left, topAligns, startState, midState);
		
		// Calculate bottom subalignments
		AlignmentList botAligns;
		align(seq1Bot, seq2Right, botAligns, midState, endState);
					
		// Combine all subalignments pairwise
		AlignmentList::const_iterator topA, botA;
		for (topA = topAligns.begin(); topA != topAligns.end(); ++topA) {
			for (botA = botAligns.begin(); botA != botAligns.end(); ++botA) {
				paths.push_back((*topA) + (*botA));
			}
		}
	}
	
	template<typename NumType>
	void
	AffineGapNWPairwiseAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2,
		  std::vector<PairwiseAlignment>& paths,
		  size_t startState,
		  size_t endState) const {
		if (seq1.size() == 0) {
			paths.push_back(PairwiseAlignment(std::string(seq2.size(), '-'),
											  seq2));
		} else if (seq2.size() == 0) {
			paths.push_back(PairwiseAlignment(seq1,
											  std::string(seq1.size(), '-')));
		} else if (seq1.size() == 1 and seq2.size() == 1) {
			align(seq1[0], seq2[0], paths, startState, endState);
		} else if (seq1.size() == 1) {
			align(seq2, seq1, paths,
				  flipState(startState), flipState(endState));
			std::for_each(paths.begin(), paths.end(),
						  std::mem_fun_ref(&PairwiseAlignment::flip));
		} else {
			// Calculate middle of first sequence
			std::string::size_type middle = seq1.size() / 2;
			std::string seq1Top = seq1.substr(0, middle);
			std::string seq1Bot = seq1.substr(middle);

			// Calculate forward values to middle of first sequence
			ScoreList forwardMatch, forwardGap1, forwardGap2;
			scorer.scoreLastRow(seq1Top, seq2,
								forwardMatch, forwardGap1, forwardGap2, startState);

			// Calculate backward values to middle of first sequence
			ScoreList backwardMatch, backwardGap1, backwardGap2;
			scorer.scoreFirstRow(seq1Bot, seq2,
								 backwardMatch, backwardGap1, backwardGap2, endState);
			
			// Combine forward and backward scores
			combineScores(forwardMatch, backwardMatch);
			combineScores(forwardGap1, backwardGap1);
			combineScores(forwardGap2, backwardGap2);

			// Find maximum positions
			typename ScoreList::iterator maxMatch, maxGap1;
			maxMatch = std::max_element(forwardMatch.begin(), forwardMatch.end());
			maxGap1 = std::max_element(forwardGap1.begin(), forwardGap1.end());

			Score maxScore = *maxMatch + *maxGap1;

			for (size_t i = 0; i < forwardMatch.size(); ++i) {
				size_t midState = ((forwardMatch[i] == maxScore ? Scorer::MATCH : 0) |
								   (forwardGap1[i] == maxScore ? Scorer::GAP1 : 0));
				if (midState != 0) {
					align(seq1, seq2, paths, startState, endState, midState,
						  middle, i);
				}
			}
		}
	}

} }

#endif // __BIO_ALIGNMENT_AFFINEGAPNWPAIRWISEALIGNER_HH__
