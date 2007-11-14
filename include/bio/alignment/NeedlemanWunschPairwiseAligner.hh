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

#ifndef __NEEDLEMAN_WUNSCH_PAIRWISE_ALIGNER_HH__
#define __NEEDLEMAN_WUNSCH_PAIRWISE_ALIGNER_HH__

#include <vector>

#include "bio/alignment/PairwiseAligner.hh"
#include "bio/alignment/NeedlemanWunschPairwiseScorer.hh"
#include "bio/alignment/ConvertingScoringMatrix.hh"
#include "math/MaxPlus.hh"
using math::MaxPlus;

namespace bio { namespace alignment {

	template<typename NumType>
	class NeedlemanWunschPairwiseAligner
		: public PairwiseAligner,
		  public NeedlemanWunschPairwiseScorer< MaxPlus<NumType> > {
	public:

		NeedlemanWunschPairwiseAligner(const ScoringMatrix<NumType>& matrix,
									   const NumType& space);
		
		PairwiseAlignment align(const std::string& seq1,
								const std::string& seq2) const;
		
		void align(const std::string& seq1,
				   const std::string& seq2,
				   AlignmentList& paths) const;
				
	private:
  		typedef typename NeedlemanWunschPairwiseScorer< MaxPlus<NumType> >::Score Score;
  		typedef typename NeedlemanWunschPairwiseScorer< MaxPlus<NumType> >::ScoreList ScoreList;

		void align(const std::string& seq1,
				   const std::string& seq2,
				   AlignmentList& paths,
				   std::string::size_type seq1Prefix,
				   std::string::size_type seq2Prefix) const;

		using NeedlemanWunschPairwiseScorer< MaxPlus<NumType> >::matrix;
		ConvertingScoringMatrix<NumType, Score> convertingMatrix;
	};

	template<typename NumType>
	NeedlemanWunschPairwiseAligner<NumType>::
	NeedlemanWunschPairwiseAligner(const ScoringMatrix<NumType>& matrix,
								   const NumType& space)
		: NeedlemanWunschPairwiseScorer< MaxPlus<NumType> >(convertingMatrix,
															space),
		  convertingMatrix(matrix) {
	}

	template<typename NumType>
	void
	NeedlemanWunschPairwiseAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2,
		  std::vector<PairwiseAlignment>& paths,
		  std::string::size_type seq1Prefix,
		  std::string::size_type seq2Prefix) const {

		std::string seq1Top = seq1.substr(0, seq1Prefix);
		std::string seq1Bot = seq1.substr(seq1Prefix);
		std::string seq2Left = seq2.substr(0, seq2Prefix);
		std::string seq2Right = seq2.substr(seq2Prefix);
		
		// Calculate top subalignments
		AlignmentList topAligns;
		align(seq1Top, seq2Left, topAligns);
		
		// Calculate bottom subalignments
		AlignmentList botAligns;
		align(seq1Bot, seq2Right, botAligns);
					
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
	NeedlemanWunschPairwiseAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2,
		  std::vector<PairwiseAlignment>& paths) const {
		if (seq1.size() == 0) {
			paths.push_back(PairwiseAlignment(std::string(seq2.size(), '-'),
											  seq2));
		} else if (seq2.size() == 0) {
			paths.push_back(PairwiseAlignment(seq1,
											  std::string(seq1.size(), '-')));
		} else if (seq1.size() == 1 and seq2.size() == 1) {
			Score spacePathValue = this->space * this->space;
			Score matchPathValue = matrix.getCharScore(seq1[0], seq2[0]);
			if (spacePathValue <= matchPathValue) {
				paths.push_back(PairwiseAlignment(seq1, seq2));
			}
			if (matchPathValue <= spacePathValue) {
				paths.push_back(PairwiseAlignment(seq1 + "-", "-" + seq2));
				paths.push_back(PairwiseAlignment("-" + seq1, seq2 + "-"));
			}
		} else if (seq1.size() == 1) {
			align(seq2, seq1, paths);
			std::for_each(paths.begin(), paths.end(),
						  std::mem_fun_ref(&PairwiseAlignment::flip));
		} else {
			// Calculate middle of first sequence
			std::string::size_type middle = seq1.size() / 2;
			std::string seq1Top = seq1.substr(0, middle);
			std::string seq1Bot = seq1.substr(middle);

			// Calculate forward values to middle of first sequence
			ScoreList forward;
			scoreLastRow(seq1Top, seq2, forward);

			// Calculate backward values to middle of first sequence
			ScoreList backward;
			std::string seq2Rev = seq2;
			std::reverse(seq2Rev.begin(), seq2Rev.end());
			std::reverse(seq1Bot.begin(), seq1Bot.end());
			scoreLastRow(seq1Bot, seq2Rev, backward);
			std::reverse(backward.begin(), backward.end());

			// Combine forward and backward scores
			std::transform(forward.begin(), forward.end(),
						   backward.begin(), forward.begin(),
						   std::multiplies<Score>());

			// Find maximum value
			Score maxElt = *std::max_element(forward.begin(), forward.end());
			
			// For each value that is maximal in middle row, calculate
			// paths that pass through that element
			for (size_t i = 0; i < forward.size(); ++i) {
				if (forward[i] == maxElt) {
					align(seq1, seq2, paths, middle, i);
				}
			}
		}
	}

	template<typename NumType>
	PairwiseAlignment
	NeedlemanWunschPairwiseAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2) const {
		if (seq1.size() == 0) {
			return PairwiseAlignment(std::string(seq2.size(), '-'), seq2);
		} else if (seq2.size() == 0) {
			return PairwiseAlignment(seq1, std::string(seq1.size(), '-'));
		} else if (seq1.size() == 1 and seq2.size() == 1) {
			Score spacePathValue = this->space * this->space;
			Score matchPathValue = matrix.getCharScore(seq1[0], seq2[0]);
			if (spacePathValue <= matchPathValue) {
				return PairwiseAlignment(seq1, seq2);
			} else {
				return PairwiseAlignment(seq1 + "-", "-" + seq2);
			}
		} else if (seq1.size() == 1) {
			PairwiseAlignment flipped = align(seq2, seq1);
			flipped.flip();
			return flipped;
		} else {
			// Calculate middle of first sequence
			std::string::size_type middle = seq1.size() / 2;
			std::string seq1Top = seq1.substr(0, middle);
			std::string seq1Bot = seq1.substr(middle);

			// Calculate forward values to middle of first sequence
			ScoreList forward;
			scoreLastRow(seq1Top, seq2, forward);

			// Calculate backward values to middle of first sequence
			ScoreList backward;
			std::string seq2Rev = seq2;
			std::reverse(seq2Rev.begin(), seq2Rev.end());
			std::reverse(seq1Bot.begin(), seq1Bot.end());
			scoreLastRow(seq1Bot, seq2Rev, backward);
			std::reverse(seq1Bot.begin(), seq1Bot.end());
			std::reverse(backward.begin(), backward.end());

			// Combine forward and backward scores
			std::transform(forward.begin(), forward.end(),
						   backward.begin(), forward.begin(),
						   std::multiplies<Score>());
			
			// Find maximum value index
			size_t seq2Prefix =
				std::max_element(forward.begin(), forward.end()) -
				forward.begin();
			std::string seq2Left = seq2.substr(0, seq2Prefix);
			std::string seq2Right = seq2.substr(seq2Prefix);
		
			// Combine best alignment of top with best alignment of bottom
			return align(seq1Top, seq2Left) +
				align(seq1Bot, seq2Right);
		}
	}
	
} }

#endif // __NEEDLEMAN_WUNSCH_PAIRWISE_ALIGNER_HH__
