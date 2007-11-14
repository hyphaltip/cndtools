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

#ifndef __CONSENSUS_ALIGNMENT_HH__
#define __CONSENSUS_ALIGNMENT_HH__

#include <string>
#include <vector>

namespace bio { namespace alignment {

	class ConsensusAlignment {
	public:
		typedef std::pair<size_t, size_t> Match;
		typedef std::vector<Match> MatchList;
		typedef MatchList::const_iterator iterator;
	
		ConsensusAlignment(const std::string& seq1, const std::string& seq2);

		void addMatch(size_t pos1, size_t pos2);
	
		ConsensusAlignment& operator+=(const ConsensusAlignment& ca);
		ConsensusAlignment& operator&=(const ConsensusAlignment& ca);
		ConsensusAlignment operator+(const ConsensusAlignment& ca) const;
		ConsensusAlignment operator&(const ConsensusAlignment& ca) const;

		const std::string& getSeq1() const;
		const std::string& getSeq2() const;

		std::pair<std::string, std::string> getAlignedSeqs() const;

		iterator begin() const;
		iterator end() const;

		void flip();
		
	private:
		std::string seq1;
		std::string seq2;
		mutable MatchList matches;
		mutable bool sorted;

		void sortMatches() const;
	
	};

	inline const std::string&
	ConsensusAlignment::getSeq1() const { return seq1; }
	inline const std::string&
	ConsensusAlignment::getSeq2() const { return seq2; }

	inline ConsensusAlignment::iterator
	ConsensusAlignment::begin() const { return matches.begin(); }
	inline ConsensusAlignment::iterator
	ConsensusAlignment::end() const { return matches.end(); }

	inline ConsensusAlignment
	ConsensusAlignment::operator+(const ConsensusAlignment& ca) const {
		ConsensusAlignment result = *this;
		return result += ca;
	}

	inline ConsensusAlignment
	ConsensusAlignment::operator&(const ConsensusAlignment& ca) const {
		ConsensusAlignment result = *this;
		return result &= ca;
	}
	
} }

#endif // __CONSENSUS_ALIGNMENT_HH__
