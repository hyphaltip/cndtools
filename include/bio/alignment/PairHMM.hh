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

#ifndef __BIO_ALIGNMENT_PAIRHMM_HH__
#define __BIO_ALIGNMENT_PAIRHMM_HH__

#include <vector>

#include "util/matrix.hh"
#include "math/LogDouble.hh"

namespace bio { namespace alignment {

    class PairHMM {
	public:
		typedef std::string Seq;
		typedef math::LogDouble Prob;
		typedef util::Matrix<Prob> ProbMatrix;
		typedef std::vector<Prob> ProbVector;
		
		PairHMM(const ProbMatrix& hEmissions,
				const ProbVector& iEmissions,
				const ProbVector& dEmissions,
				const ProbVector& beginProbs,
				const ProbVector& endProbs,
				const ProbMatrix& transitions);

		Prob likelihood(const Seq& s1,
						const Seq& s2);

		enum State { H_STATE, I_STATE, D_STATE };
		
	private:
		const ProbMatrix hEmissions;
		const ProbVector iEmissions;
		const ProbVector dEmissions;
		const ProbVector& beginProbs;
		const ProbVector& endProbs;
		const ProbMatrix transitions;
    };

} }

#endif // __BIO_ALIGNMENT_PAIRHMM_HH__
