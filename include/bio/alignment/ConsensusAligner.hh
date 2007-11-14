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

#ifndef __CONSENSUS_ALIGNER_HH__
#define __CONSENSUS_ALIGNER_HH__

#include <vector>

#include "bio/alignment/ConsensusAlignment.hh"
#include "bio/alignment/NeedlemanWunschPairwiseScorer.hh"
#include "math/MaxPlus.hh"
using math::MaxPlus;

namespace bio { namespace alignment {

	class ConsensusAligner<typename Score> :
		private NeedlemanWunschPairwiseScorer< MaxPlus<Score> > {
	public:
		ConsensusAligner(const Score& match,
						 const Score& mismatch,
						 const Score& space);
		
		ConsensusAlignment align(const std::string& seq1,
								 const std::string& seq2) const;
	
	protected:
		ConsensusAlignment align(const std::string& seq1,
								 const std::string& seq2,
								 std::string::size_type seq1Prefix,
								 std::string::size_type seq2Prefix);
	};

	template<typename Score>
	ConsensusAligner<Score>::
	ConsensusAligner(const Score& match,
					 const Score& mismatch,
					 const Score& space)
		: NeedlemanWunschPairwiseScorer< MaxPlus<Score> >(match, mismatch, space)
	{}
	
	template<typename Score>
	ConsensusAlignment
	ConsensusAligner<Score>::
	align(const std::string& seq1,
		  const std::string& seq2,
		  std::string::size_type seq1Prefix,
		  std::string::size_type seq2Prefix) const {
		std::string seq1Top = seq1.substr(0, seq1Prefix);
		std::string seq1Bot = seq1.substr(seq1Prefix);
		std::string seq2Left = seq2.substr(0, seq2Prefix);
		std::string seq2Right = seq2.substr(seq2Prefix);

		return align(seq1Top, seq2Left) + align(seq1Bot, seq2Right);
	}
		
	template<typename Score>
	ConsensusAlignment
	align(const std::string& seq1,
		  const std::string& seq2) const {
		if (seq1.size() == 0 or seq2.size() == 0) {
			return ConsensusAlignment(seq1, seq2);
		} else if (seq1.size() == 1 and seq2.size() == 1) {
			Score spacePathValue = space * space;
			Score matchPathValue = (seq1[0] == seq2[0] ? match : mismatch);
			ConsensusAlignment ca(seq1, seq2);			
			if (spacePathValue < matchPathValue) {
				ca.addMatch(0, 0);
			}
			return ca;
		} else if (seq1.size() == 1) {
			ConsensusAlignment ca = alignGlobal(seq2, seq1);
			ca.flip();
			return ca;
		} else {
			// Calculate middle of first sequence
			std::string::size_type middle = seq1.size() / 2;
			std::string seq1Top = seq1.substr(0, middle);
			std::string seq1Bot = seq1.substr(middle);

			// Calculate forward values to middle of first sequence
			std::vector<Score> forward;
			pairwiseAligner.alignGlobalLastRow(seq1Top, seq2, forward);

			// Calculate backward values to middle of first sequence
			std::vector<Score> backward;
			std::string seq2Rev = seq2;
			std::reverse(seq2Rev.begin(), seq2Rev.end());
			std::reverse(seq1Bot.begin(), seq1Bot.end());
			pairwiseAligner.alignGlobalLastRow(seq1Bot, seq2Rev, backward);
			std::reverse(backward.begin(), backward.end());

			// Combine forward and backward scores
			std::transform(forward.begin(), forward.end(),
						   backward.begin(), forward.begin(),
						   std::multiplies<Score>());

			// Find maximum value
			Score maxElt = *std::max_element(forward.begin(), forward.end());

			ConsensusAlignment ca(seq1, seq2);
			bool firstAlignment = true;
			// For each value that is maximal in middle row, calculate
			// paths that pass through that element
			for (size_t i = 0; i < forward.size(); ++i) {
				if (forward[i] == maxElt) {
					if (firstAlignment) {
						ca = alignGlobal(seq1, seq2, middle, i);
					} else {
						ca &= alignGlobal(seq1, seq2, middle, i);
					}
				}
			}

			return ca;
		}
	}

} }


#endif // __CONSENSUS_ALIGNER_HH__
