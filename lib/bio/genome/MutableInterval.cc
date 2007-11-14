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

#include <stdexcept>

#include "bio/genome/MutableInterval.hh"
#include "util/string.hh"

namespace bio { namespace genome {

	void MutableInterval::setStartCoord(const Coord& c) {
		setChrom(c.chrom);
		setStart(c.pos);
	}
	
	void MutableInterval::setEndCoord(const Coord& c) {
		setChrom(c.chrom);
		setEnd(c.pos);
	}

	void MutableInterval::setInterval(const Interval& i) {
		setChrom(i.getChrom());
		setStart(i.getStart());
		setEnd(i.getEnd());
		setStrand(i.getStrand());
	}

	MutableInterval& MutableInterval::operator=(const Interval& i) {
		setInterval(i);
		return *this;
	}

	void MutableInterval::flip() { setStrand(getStrand().opposite()); }

	std::istream& operator>>(std::istream& strm, MutableInterval& i) {
		std::string s;
		strm >> s;
		if (not strm) { return strm; }
		size_t colonPos = s.find(':');
		if (colonPos == std::string::npos) {
			throw std::runtime_error("Invalid interval: " + s);
		}
		size_t dashPos = s.find('-', colonPos);
		if (dashPos == std::string::npos) {
			throw std::runtime_error("Invalid interval: " + s);
		}
		size_t strandPos = s.size() - 1;
		if (strandPos <= dashPos) {
			throw std::runtime_error("Invalid interval: " + s);
		}
		util::string::Converter<genome::Position> toPosition;
		i.setChrom(s.substr(0, colonPos));
		i.setStart(toPosition(s.substr(colonPos + 1, dashPos - colonPos - 1)));
		i.setEnd(toPosition(s.substr(dashPos + 1, strandPos - dashPos - 1)));
		i.setStrand(s[strandPos]);
		return strm;
	}
	
} }
