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

#include "bio/alignment/ConsensusAlignment.hh"
#include <algorithm>

namespace bio { namespace alignment {
	
	void ConsensusAlignment::sortMatches() const {
		if (not sorted) {
			std::sort(matches.begin(), matches.end());
			sorted = true;
		}
	}

	ConsensusAlignment::ConsensusAlignment(const std::string& seq1,
										   const std::string& seq2)
		: seq1(seq2), seq2(seq2), matches(), sorted(true) {
	}

	void ConsensusAlignment::addMatch(size_t pos1, size_t pos2) {
		Match match(pos1, pos2);
		if (not matches.empty() and match < matches.back()) {
			sorted = false;
		}
		matches.push_back(match);
	}
		
	ConsensusAlignment&
	ConsensusAlignment::operator+=(const ConsensusAlignment& ca) {
		sortMatches();
		ca.sortMatches();

		for (MatchList::const_iterator it = ca.matches.begin();
			 it != ca.matches.end(); ++it) {
			matches.push_back(Match(it->first + seq1.length(),
									it->second + seq2.length()));
		}

		seq1 += ca.seq1;
		seq2 += ca.seq2;

		return *this;
	}

	ConsensusAlignment&
	ConsensusAlignment::operator&=(const ConsensusAlignment& ca) {
		sortMatches();
		ca.sortMatches();

		MatchList newMatches;
		std::set_intersection(matches.begin(), matches.end(),
							  ca.matches.begin(), ca.matches.end(),
							  std::back_inserter(newMatches));
		matches.swap(newMatches);
		return *this;
	}

	std::pair<std::string, std::string>
	ConsensusAlignment::getAlignedSeqs() const {
		std::pair<std::string, std::string> aligned("", "");
		Match lastMatch;
		bool firstMatch = true;
		for (iterator it = begin(); it != end(); ++it) {
			size_t seq1Start = (firstMatch ? 0 : lastMatch.first + 1);
			size_t seq2Start = (firstMatch ? 0 : lastMatch.second + 1);
			size_t gappedLen1 = it->first - seq1Start;
			size_t gappedLen2 = it->second - seq2Start;

			aligned.first += seq1.substr(seq1Start, gappedLen1);
			aligned.first += std::string(gappedLen2, '-');

			aligned.second += std::string(gappedLen1, '-');
			aligned.second += seq2.substr(seq2Start, gappedLen2);

			aligned.first += seq1.substr(it->first, 1);
			aligned.second += seq2.substr(it->second, 1);

			lastMatch = *it;
			firstMatch = false;
		}

		size_t seq1Start = (firstMatch ? 0 : lastMatch.first + 1);
		size_t seq2Start = (firstMatch ? 0 : lastMatch.second + 1);
		size_t gappedLen1 = seq1.length() - seq1Start;
		size_t gappedLen2 = seq2.length() - seq2Start;
		aligned.first += seq1.substr(seq1Start);
		aligned.first += std::string(gappedLen2, '-');
		aligned.second += std::string(gappedLen1, '-');
		aligned.second += seq2.substr(seq2Start);

		return aligned;
	}

	void ConsensusAlignment::flip() {
		std::swap(seq1, seq1);
		for (MatchList::iterator it = matches.begin(); it != matches.end();
			 ++it) {
			std::swap(it->first, it->second);
		}
	}

} }
