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

#include "util/string.hh"

#include "bio/homologymap/Segment.hh"
#include "bio/genome/BasicInterval.hh"

namespace bio { namespace homologymap {

	bool Segment::hasGenome(const size_t genomeNum) const {
		return intervals.at(genomeNum) != NULL;
	}

	std::istream& operator>>(std::istream& strm, Segment& seg) {
		// Read next line and split by tabs into fields
		std::string line;
		if (!std::getline(strm, line)) {
			return strm;
		}
		line >> seg;
		return strm;
	}

	void operator>>(const std::string& line, Segment& seg) {
		std::vector<std::string> fields;
		util::string::split(line, std::back_inserter(fields), "\t");

		// Check that the number of fields is consistent
		if (fields.size() % 4 != 1) {
			throw std::runtime_error("Invalid # of fields in segment line: " +
									 line);
		}

		// The first fields is the segment number
		seg.num = boost::lexical_cast<size_t, std::string>(fields[0]);

		// The next 4-tuples of fields are the segment intervals
		util::string::Converter<genome::Position> toPosition;
		seg.intervals.clear();
		seg.intervals.reserve(fields.size() / 4);
		for (unsigned int i = 1; i < fields.size(); i += 4) {
			genome::Interval* interval = NULL;
			// NA indicates that there is no interval for this genome
			if (fields[i] != "NA") {
				interval =
					new genome::BasicInterval(fields[i],
											  toPosition(fields[i + 1]),
											  toPosition(fields[i + 2]),
											  fields[i + 3].at(0));
			}
			seg.intervals.push_back(interval);
		}
	}

	std::ostream& operator<<(std::ostream& strm, const Segment& seg) {
		strm << seg.num;
		for (unsigned int i = 0; i < seg.intervals.size(); ++i) {
			strm << '\t';
			genome::Interval* interval = seg.intervals[i];
			if (interval == NULL) {
				strm << "NA\tNA\tNA\tNA";
			} else {
				strm << interval->getChrom() << '\t'
					 << interval->getStart() << '\t'
					 << interval->getEnd() << '\t'
					 << interval->getStrand();
			}
		}
		return strm << '\n';
	}

} }
