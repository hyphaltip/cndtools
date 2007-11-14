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

#ifndef __BIO_ALIGNMENT_ALIGNMENTSLICER_HH__
#define __BIO_ALIGNMENT_ALIGNMENTSLICER_HH__

#include <string>
#include <vector>
#include <map>

#include "filesystem.hh"
#include "bio/alignment/BasicNamedMultipleAlignment.hh"
#include "bio/genome/Interval.hh"
#include "bio/homologymap/Map.hh"

namespace bio { namespace alignment {

	class AlignmentSlicer {
	public:
		AlignmentSlicer(const filesystem::Path& alignDir,
						const std::string& sourceGenome,
						char noMapChar = '-');
	
		BasicNamedMultipleAlignment
		getSlice(const bio::genome::Interval& sourceInterval);

	private:
		bool readAlignment(size_t segNum);

		static const std::string GENOMES_FILENAME;
		static const std::string MAP_FILENAME;
		
		filesystem::Path alignDir;
		std::string sourceGenome;
		char noMapChar;
		std::vector<std::string> genomes;
		std::map<std::string, size_t> indexMap;
		size_t sourceIndex;
		homologymap::Map map;
		BasicNamedMultipleAlignment align;
		size_t lastSegNum;
	};

} }

#endif // __BIO_ALIGNMENT_ALIGNMENTSLICER_HH__
