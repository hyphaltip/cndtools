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

#include "bio/formats/agp/Record.hh"

namespace bio { namespace formats { namespace agp {

	// Source sequence constructor
	Record::Record(const std::string& chrom,
				   const std::string& contig,
				   const genome::Position contigStart,
				   const genome::Position contigEnd,
				   const unsigned int agpRecNum,
				   const char sourceType,
				   const std::string& sourceAccession,
				   const genome::Position sourceStart,
				   const genome::Position sourceEnd,
				   const char orientation)
		: fields(9, "") {
		if (sourceType == 'N') {
			throw std::runtime_error("Attempted construction of non-gap AGP "
									 "record with source type = N");
		}
		fields[0] = (contig.empty() ? chrom : chrom + "/" + contig);
		fields[1] = util::string::toString(contigStart);
		fields[2] = util::string::toString(contigEnd);
		fields[3] = util::string::toString(agpRecNum);
		fields[4] = util::string::toString(sourceType);
		fields[5] = sourceAccession;
		fields[6] = util::string::toString(sourceStart);
		fields[7] = util::string::toString(sourceEnd);
		fields[8] = util::string::toString(orientation);
	}
	
	// Gap constructor
	Record::Record(const std::string& chrom,
				   const std::string& contig,
				   const genome::Position contigStart,
				   const genome::Position contigEnd,
				   const unsigned int agpRecNum,						  
				   const char sourceType,
				   const genome::Position gapLength,
				   const std::string& gapKind,
				   const bool gapBridged)
		: fields(8, "") {
		if (sourceType != 'N') {
			throw std::runtime_error("Attempted construction of gap AGP "
									 "record with source type != N");
		}
		fields[0] = (contig.empty() ? chrom : chrom + "/" + contig);
		fields[1] = util::string::toString(contigStart);
		fields[2] = util::string::toString(contigEnd);
		fields[3] = util::string::toString(agpRecNum);
		fields[4] = util::string::toString(sourceType);
		fields[5] = util::string::toString(gapLength);
		fields[6] = gapKind;
		fields[7] = (gapBridged ? "yes" : "no");
	}		

	void Record::stripComment(std::string& line) {
		std::string::size_type pos = line.find('#');
		if (pos == std::string::npos) {
			return;
		} else {
			while (pos > 0 and std::isspace(line[pos - 1])) {
				--pos;
			}
			line.erase(pos);
		}
	}
	
	std::ostream& operator<<(std::ostream& strm, const Record& r) {
		std::vector<std::string>::const_iterator it = r.fields.begin();
		strm << *it;
		for (++it; it != r.fields.end(); ++it) {
			strm << '\t' << *it;
		}
		return strm << '\n';
	}
	
	std::istream& operator>>(std::istream& strm, Record& r) {
		std::string line;
		while (line.empty() and std::getline(strm, line)) {
			Record::stripComment(line);
		}
		if (not line.empty()) {
			line >> r;
		}
		return strm;
	}

	void operator>>(const std::string& line, Record& r) {
		r.fields.clear();
		util::string::split(line, std::back_inserter(r.fields), "\t");
		if (not ((r.fields.size() == 8 and r.isGap()) ||
				 (r.fields.size() == 9 and not r.isGap()))) {
			throw std::runtime_error("Invalid AGP record line:\n" + line);
		}
	}

} } }
