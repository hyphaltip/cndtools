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

#include "boost/lexical_cast.hpp"

#include "util/string.hh"
#include "bio/repeatmasker/RepeatMaskerRecord.hh"

namespace bio { namespace repeatmasker {

	RepeatMaskerRecord::RepeatMaskerRecord() {
		fields.reserve(16);
	}
	
	int RepeatMaskerRecord::getScore() const {
		return boost::lexical_cast<int>(fields[0]);
	}
		
	float RepeatMaskerRecord::getPctDivergence() const {
		return boost::lexical_cast<float>(fields[1]);		
	}
	
	float RepeatMaskerRecord::getPctDeleted() const {
		return boost::lexical_cast<float>(fields[2]);
	}
	
	float RepeatMaskerRecord::getPctInserted() const {
		return boost::lexical_cast<float>(fields[3]);
	}

	genome::BasicInterval RepeatMaskerRecord::getQueryInterval() const {
		return genome::BasicInterval(getQueryName(),
									 getQueryStart(),
									 getQueryEnd());
	}
	
	std::string RepeatMaskerRecord::getQueryName() const {
		return fields[4];
	}
	
	genome::Position RepeatMaskerRecord::getQueryStart() const {
		return boost::lexical_cast<genome::Position>(fields[5]);
	}
	
	genome::Position RepeatMaskerRecord::getQueryEnd() const {
		return boost::lexical_cast<genome::Position>(fields[6]);
	}
	
	genome::Distance RepeatMaskerRecord::getQueryLeft() const {
		return boost::lexical_cast<genome::Distance>(fields[7]);
	}
	
	bool RepeatMaskerRecord::isComplement() const {
		return fields[8][0] == 'C';
	}
	
	std::string RepeatMaskerRecord::getRepeatClass() const {
		return fields[10];
	}

	genome::BasicInterval RepeatMaskerRecord::getRepeatInterval() const {
		return genome::BasicInterval(getRepeatName(),
									 getRepeatStart(),
									 getRepeatEnd(),
									 (isComplement() ? '-' : '+'));
	}
	
	std::string RepeatMaskerRecord::getRepeatName() const {
		return fields[9];
	}
	
	genome::Position RepeatMaskerRecord::getRepeatStart() const {
		if (isComplement()) {
			return boost::lexical_cast<genome::Position>(fields[13]);
		} else {
			return boost::lexical_cast<genome::Position>(fields[11]);
		}
	}
	
	genome::Position RepeatMaskerRecord::getRepeatEnd() const {
		return boost::lexical_cast<genome::Position>(fields[12]);
	}
	
	genome::Distance RepeatMaskerRecord::getRepeatLeft() const {
		if (isComplement()) {
			return boost::lexical_cast<genome::Distance>(fields[11]);
		} else {
			return boost::lexical_cast<genome::Distance>(fields[13]);
		}
	}
	
	int RepeatMaskerRecord::getID() const {
		return boost::lexical_cast<int>(fields[14]);
	}
	
	bool RepeatMaskerRecord::isIncludedInHigherScoringMatch() const {
		return fields.size() == 16;
	}		
	
	std::ostream& operator<<(std::ostream& strm, const RepeatMaskerRecord& r) {
		std::vector<std::string>::const_iterator it = r.fields.begin();
		strm << *it;
		for (++it; it != r.fields.end(); ++it) {
			strm << '\t' << *it;
		}
		return strm << '\n';
	}
	
	std::istream& operator>>(std::istream& strm, RepeatMaskerRecord& r) {
		std::string line;
		while (line.empty() and std::getline(strm, line)) {}
		if (not line.empty()) {
			line >> r;
		}
		return strm;
	}

	void operator>>(const std::string& line, RepeatMaskerRecord& r) {
		r.fields.clear();
		util::string::split(line, std::back_inserter(r.fields));
		if (not (r.fields.size() == 15 or r.fields.size() == 16)) {
			throw std::runtime_error("Invalid RepeatMasker record line:\n" + line);
		}
	}

} }
