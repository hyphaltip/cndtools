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

#include "bio/alignment/PairHMM.hh"

namespace bio { namespace alignment {

		PairHMM::PairHMM(const ProbMatrix& hEmissions,
						 const ProbVector& iEmissions,
						 const ProbVector& dEmissions,
						 const ProbVector& beginProbs,
						 const ProbVector& endProbs,
						 const ProbMatrix& transitions) :
			hEmissions(hEmissions),
			iEmissions(iEmissions),
			dEmissions(dEmissions),
			beginProbs(beginProbs),
			endProbs(endProbs),
			transitions(transitions) {
		}

		PairHMM::Prob
		PairHMM::likelihood(const Seq& seq1,
							const Seq& seq2) {
			const Prob zero = Prob(0);
			const Prob one = Prob(1);

			// Initialize vectors
			ProbVector h_state(seq2.size() + 1, zero);
			ProbVector i_state(seq2.size() + 1, zero);
			ProbVector d_state(seq2.size() + 1, zero);

			// Calculate first row
			if (seq2.size() > 0) {
				i_state[1] = beginProbs[I_STATE];
				i_state[1] *= iEmissions[seq2[0]];
				for (size_t j = 2; j <= seq2.size(); ++j) {
					i_state[j] = transitions(I_STATE, I_STATE) * i_state[j - 1];
					i_state[j] *= iEmissions[seq2[j - 1]];
				}
			}
			
			// Calculate all rows
			Prob up_h_state, up_d_state, up_i_state;
			Prob diag_h_state, diag_d_state, diag_i_state;
			for (size_t i = 1; i <= seq1.size(); ++i) {
				for (size_t j = 0; j <= seq2.size(); ++j) {
					// Keep entries from last row
					up_h_state = h_state[j];
					up_i_state = i_state[j];
					up_d_state = d_state[j];

					// Calculate H and I states
					if (j == 0) {
						h_state[0] = i_state[0] = zero;
					} else {
						h_state[j] = (i == 1 and j == 1) ?
							beginProbs[H_STATE] :
							transitions(H_STATE, H_STATE) * diag_h_state +
							transitions(I_STATE, H_STATE) * diag_i_state +
							transitions(D_STATE, H_STATE) * diag_d_state;
						h_state[j] *= hEmissions(seq1[i - 1], seq2[j - 1]);

						i_state[j] =
							transitions(H_STATE, I_STATE) * h_state[j - 1] +
							transitions(I_STATE, I_STATE) * i_state[j - 1] +
							transitions(D_STATE, I_STATE) * d_state[j - 1];
						i_state[j] *= iEmissions[seq2[j - 1]];
					}

					// Calculate D state
					d_state[j] = (i == 1 and j == 0) ?
						beginProbs[D_STATE] :
						transitions(H_STATE, D_STATE) * up_h_state +
						transitions(I_STATE, D_STATE) * up_i_state +
						transitions(D_STATE, D_STATE) * up_d_state;
					d_state[j] *= dEmissions[seq1[i - 1]];

					// Entries from last row become diagonal entries
					diag_h_state = up_h_state;
					diag_i_state = up_i_state;
					diag_d_state = up_d_state;
				}

			}

			return
				h_state[seq2.size()] * endProbs[H_STATE] +
				i_state[seq2.size()] * endProbs[I_STATE] +
				d_state[seq2.size()] * endProbs[D_STATE];
		}
	
	} }
