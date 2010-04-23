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

#ifndef __BIO_HOMOLOGYMAP_HOMOLOGYMAPPER_HH__
#define __BIO_HOMOLOGYMAP_HOMOLOGYMAPPER_HH__

#include <string>
#include <vector>

#include "bio/genome/IntervalMapper.hh"
#include "boost/unordered_map.hpp"
#include "bio/genome/Interval.hh"
#include "bio/alignment/BasicNamedMultipleAlignment.hh"
#include "bio/homologymap/Map.hh"
#include "filesystem.hh"

namespace bio { namespace homologymap {

    class HomologyMapper : public genome::IntervalMapper {
	public:
		
		HomologyMapper(const filesystem::Path& alignDir,
					   const std::string& sourceGenome,
					   const std::string& targetGenome);

		void map(const genome::Interval& i,
				 std::vector<genome::BasicInterval>& mapped);
		
		void map(const genome::Interval& i,
				 std::vector<genome::BasicInterval>& mapped,
				 const std::string& source,
				 const std::string& target);

		void setTargetGenome(const std::string& targetGenome);
		void setSourceGenome(const std::string& sourceGenome);

	private:
		size_t getIndex(const std::string& g) const;
		bool readAlignment(size_t segmentNum);	
		void flipInterval(alignment::Interval& i, size_t length);
		void getInterval(genome::Interval* segSourceInt,
						 genome::Interval* segTargetInt,
						 int alignSourceNum,
						 int alignTargetNum,
						 const genome::Interval& sourceInterval,
						 genome::BasicInterval& targetInterval);

		static const std::string GENOMES_FILENAME;
		static const std::string MAP_FILENAME;
		
		Map hmap;
		filesystem::Path alignDir;
		alignment::BasicNamedMultipleAlignment align;
		size_t lastSegmentNum;
        boost::unordered_map<std::string, size_t> genomeIndices;
		std::string sourceGenome;
		std::string targetGenome;
	};
		
} }

#endif // __BIO_HOMOLOGYMAP_HOMOLOGYMAPPER_HH__
