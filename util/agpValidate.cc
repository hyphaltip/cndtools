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

#include <iostream>
#include <stdexcept>
#include <vector>
#include <map>

#include "bio/formats/agp/Record.hh"
#include "bio/formats/agp/InputStream.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "util/stl.hh"
using namespace bio::formats::agp;
using util::string::toString;
using util::stl::hash_set;
using util::stl::hash_map;


class AGPValidationError : public std::runtime_error {
public:
	AGPValidationError(const std::string& what, const Record& rec) :
		std::runtime_error(what + " in record:\n" + toString(rec)) {
	}
};

struct AGPContigStartSorter {
	bool operator()(const Record* r1, const Record* r2) const {
		return r1->getContigStart() < r2->getContigStart();
	}
};

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	std::string description =
		"Validate AGP input";
	util::options::Parser parser("< agpInput", description);
	parser.parse(argv, argv + argc);

	try {
		typedef std::vector<Record*> RecList;
		typedef std::map<std::string, RecList> ChromRecsMap;
		ChromRecsMap chromRecs;
		
		// Read AGP records into map
		InputStream agpStream(std::cin);
		Record rec;
		while (agpStream >> rec) {
			chromRecs[rec.getChrom()].push_back(new Record(rec));
		}
		
		ChromRecsMap::iterator it;
		for (it = chromRecs.begin(); it != chromRecs.end(); ++it) {
			const std::string& chrom = it->first;
			RecList& recs = it->second;
			
			// Order the AGP records by contig start coordinate
			std::sort(recs.begin(), recs.end(), AGPContigStartSorter());
			
			bio::genome::Position lastEnd = 0;
			RecList::const_iterator recIt;
			for (recIt = recs.begin(); recIt != recs.end(); ++recIt) {
				Record& currRec = **recIt;

				// Check that contig coordinates are valid
				if (currRec.getContigStart() > currRec.getContigEnd()) {
					throw AGPValidationError("Contig start is greater than "
											 "contig end", currRec);
				}

				// Check that length of source is the same as contig length
				if (currRec.isGap()) {
					if (currRec.getGapLength() != currRec.getContigLength()) {
						throw AGPValidationError("Gap length is not equal to "
												 "contig length", currRec);
					}
				} else {
					// Check that source coordinates are valid
					if (currRec.getSourceStart() > currRec.getSourceEnd()) {
						throw AGPValidationError("Source start is greater "
												 "than source end", currRec);
					}
					if (currRec.getSourceLength() != currRec.getContigLength()) {
						throw AGPValidationError("Source length is not equal "
												 "to contig length", currRec);
					}
				}

				// Check that the start in the contig is one more than
				// the last end
				if (currRec.getContigStart() != lastEnd + 1) {
					throw AGPValidationError("Gap or overlap in contig",
											 currRec);
				}

				lastEnd = currRec.getContigEnd();
			}

			std::cout << chrom << '\t' << lastEnd << '\n';
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
