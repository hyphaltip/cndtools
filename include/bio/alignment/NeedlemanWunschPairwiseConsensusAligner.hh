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

#ifndef __BIO_ALIGNMENT_NEEDLEMANWUNSCHPAIRWISECONSENSUSALIGNER_HH__
#define __BIO_ALIGNMENT_NEEDLEMANWUNSCHPAIRWISECONSENSUSALIGNER_HH__

#include <vector>

#include "bio/alignment/PairwiseAligner.hh"
#include "bio/alignment/NeedlemanWunschPairwiseScorer.hh"
#include "bio/alignment/ConvertingScoringMatrix.hh"
#include "math/MaxPlus.hh"
#include "util/string.hh"
using math::MaxPlus;

namespace bio { namespace alignment {

	template<typename NumType>
	class NeedlemanWunschPairwiseConsensusAligner
		: public PairwiseAligner,
		  protected NeedlemanWunschPairwiseScorer< MaxPlus<NumType> > {
	public:
		
		NeedlemanWunschPairwiseConsensusAligner(const ScoringMatrix<NumType>& matrix,
												const NumType& space);
		
		PairwiseAlignment align(const std::string& seq1,
								const std::string& seq2) const;
				
	protected:
		typedef typename NeedlemanWunschPairwiseScorer< MaxPlus<NumType> >::Score Score;
		typedef typename NeedlemanWunschPairwiseScorer< MaxPlus<NumType> >::ScoreMatrix ScoreMatrix;

		void findConsensusEdges(const std::string& seq1,
								const std::string& seq2,
								const ScoreMatrix& forward,
								const ScoreMatrix& backward,
								std::vector<size_t>& alignment,
								std::vector<bool>& isAligned) const;

		PairwiseAlignment constructAlignment(const std::string& seq1,
											 const std::string& seq2,
											 const std::vector<size_t>& alignment1,
											 const std::vector<bool>& isAligned1,
											 const std::vector<size_t>& alignment2,
											 const std::vector<bool>& isAligned2) const;

		ConvertingScoringMatrix<NumType, Score> convertingMatrix;
	};

	template<typename NumType>
	NeedlemanWunschPairwiseConsensusAligner<NumType>::
	NeedlemanWunschPairwiseConsensusAligner(const ScoringMatrix<NumType>& matrix,
											const NumType& space)
		: NeedlemanWunschPairwiseScorer< MaxPlus<NumType> >(convertingMatrix,
															space),
		  convertingMatrix(matrix) {
	}

	template<typename NumType>
	void
	NeedlemanWunschPairwiseConsensusAligner<NumType>::
	findConsensusEdges(const std::string& seq1,
					   const std::string& seq2,
					   const ScoreMatrix& forward,
					   const ScoreMatrix& backward,
					   std::vector<size_t>& alignment,
					   std::vector<bool>& isAligned) const {
		// Initialize alignment vectors
		alignment.resize(forward.getNumRows() - 1);
		isAligned.resize(forward.getNumRows() - 1);		
		std::fill(alignment.begin(), alignment.end(), forward.getNumCols() - 1);
		std::fill(isAligned.begin(), isAligned.end(), false);
		
		Score optimal = backward(0, 0);
			
		for (size_t i = 0; i < forward.getNumRows() - 1; ++i) {
			for (size_t j = 0; j < forward.getNumCols(); ++j) {
				if ((forward(i, j) * backward(i, j)) == optimal) {
					bool isGapPathOpt = (backward(i, j) ==
										 space * backward(i + 1, j));
					bool isMatchPathOpt = j < (forward.getNumCols() - 1) and
						(backward(i, j) ==
						 matrix.getCharScore(seq1[i], seq2[j]) * backward(i + 1, j + 1));
					
					if ((isGapPathOpt and isMatchPathOpt) or
						(isMatchPathOpt and isAligned[i]) or
						(isGapPathOpt and alignment[i] != m)) {
						isAligned[i] = false;
						break;
					} else if (isGapPathOpt or isMatchPathOpt) {
						isAligned[i] = true;
						if (isMatchPathOpt) {
							alignment[i] = j;
						}
					}
				}
			}
		}
	}

	template<typename NumType>
	PairwiseAlignment
	NeedlemanWunschPairwiseConsensusAligner<NumType>::
	constructAlignment(const std::string& seq1,
					   const std::string& seq2,
					   const std::vector<size_t>& alignment1,
					   const std::vector<bool>& isAligned1,
					   const std::vector<size_t>& alignment2,
					   const std::vector<bool>& isAligned2) const {		
		std::string alignString1, alignString2;
		size_t i = 0, j = 0;
		while (i < seq1.size() or j < seq2.size()) {
			if (i != seq1.size() and
				(j == seq2.size()
				 or (not isAligned1[i])
				 or alignment1[i] == seq2.size())) {
				alignString1.push_back(seq1.at(i));
				alignString2.push_back(isAligned1[i] ? '-' : '=');
				++i;
			} else if (i == seq1.size() or (not isAligned2[j])
					   or alignment2[j] == seq1.size()) {
				alignString1.push_back(isAligned2[j] ? '-' : '=');
				alignString2.push_back(seq2.at(j));
				++j;
			} else {
				alignString1.push_back(seq1.at(i));
				alignString2.push_back(seq2.at(j));
				++i;
				++j;
			}
		}
		
		return PairwiseAlignment(alignString1, alignString2);
	}
	
	template<typename NumType>
	PairwiseAlignment
	NeedlemanWunschPairwiseConsensusAligner<NumType>::
	align(const std::string& seq1,
		  const std::string& seq2) const {
		using util::string::reverse_copy;
			
		ScoreMatrix forward, backward;

		scoreMatrix(seq1, seq2, forward);
		scoreMatrix(reverse_copy(seq1), reverse_copy(seq2), backward);
		backward.reverse();

		std::vector<bool> isAligned1, isAligned2;
		std::vector<size_t> alignment1, alignment2;

		findConsensusEdges(seq1, seq2, forward, backward,
						   alignment1, isAligned1);
		forward.transpose();
		backward.transpose();
		findConsensusEdges(seq2, seq1,
						   forward, backward, alignment2, isAligned2);

		return constructAlignment(seq1, seq2,
								  alignment1, isAligned1,
								  alignment2, isAligned2);
	}
	
} }

#endif // __BIO_ALIGNMENT_NEEDLEMANWUNSCHPAIRWISECONSENSUSALIGNER_HH__
