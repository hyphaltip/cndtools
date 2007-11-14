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

#ifndef __BIO_GENOME_AGPFORWARDMAPPER_HH__
#define __BIO_GENOME_AGPFORWARDMAPPER_HH__

#include <vector>
#include <list>
#include <istream>

#include "bio/formats/agp/Record.hh"
#include "bio/genome/IntervalMapper.hh"
#include "util/stl.hh"

namespace bio { namespace agp {

	class AGPForwardMapper : public genome::IntervalMapper {
	public:
		AGPForwardMapper(std::istream& stream);
		void map(const genome::Interval& i,
				 std::vector<genome::BasicInterval>& mapped);
		genome::Distance getChromSize(const std::string& chrom) const;

	private:
		bool overlap(const genome::Interval& i,
					 const formats::agp::Record& rec);
		
		typedef std::list<formats::agp::Record> AGPList;
		typedef util::stl::hash_map<std::string, AGPList> AGPMap;
		typedef util::stl::hash_map<std::string, genome::Distance> LengthMap;
		AGPMap contigMap;
		LengthMap chromSizes;
	};

} }

#endif // __BIO_GENOME_AGPFORWARDMAPPER_HH__
