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
#include <vector>
#include <map>

#include "bio/formats/agp.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "util/stl.hh"
using namespace bio::formats::agp;
using util::string::toString;
using util::stl::hash_set;
using util::stl::hash_map;

struct AGPSorter {
	bool operator()(const Record2* r1, const Record2* r2) const {
		return *r1 < *r2;
	}
};

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	std::vector<std::string> inconsistentChromsVect;
	
	std::string description =
		"Creates an AGP mapping from the objects in the input AGP file to "
		"objects that are consistent (have components that are ordered and "
		"oriented)";
	util::options::Parser parser("< agpInput", description);
	parser.addAppendArg("inconsistentChrom",
						"Name of object (chromosome) in input that is "
						"inconsistent", inconsistentChromsVect);
	parser.parse(argv, argv + argc);

	try {
		hash_set<std::string> inconsistentChroms(inconsistentChromsVect.begin(),
												 inconsistentChromsVect.end());

		typedef std::vector<Record2*> RecList;
		typedef std::map<std::string, RecList> ChromRecsMap;
		ChromRecsMap chromRecs;
	
		// Read AGP records into map
		InputStream agpInStream(std::cin);
		OutputStream agpOutStream(std::cout);
		Record2 rec;
		while (agpInStream >> rec) {
			chromRecs[rec.getChrom()].push_back(new Record2(rec));
		}

		size_t recNum = 1;
		ChromRecsMap::iterator it;
		for (it = chromRecs.begin(); it != chromRecs.end(); ++it) {
			const std::string& chrom = it->first;
			RecList& recs = it->second;

			// Order the AGP records by contig start coordinate
			std::sort(recs.begin(), recs.end(), AGPSorter());

			if (inconsistentChroms.find(chrom) == inconsistentChroms.end()) {
				// Chromosome is consisent.  Output a single record
				// spanning entire chromosome
				bio::genome::Distance len = recs.back()->getEnd();
				Record2 rec;
				rec.setChrom(chrom);
				rec.setStart(0);
				rec.setEnd(len);
				rec.partNumber = recNum;
				rec.componentType = Constants::OTHER;
				rec.setSourceChrom(chrom);
				rec.setSourceStart(0);
				rec.setSourceEnd(len);
				rec.setSourceStrand('+');
				agpOutStream << rec;
				++recNum;
			} else {
				// Chromosome is inconsistent.  Output consistent source
				// records (contigs with spanned gaps) as single records
				size_t chromRecNum = 1;
				Record2 lastRec;
				lastRec.setChrom("");
				lastRec.setStart(0);
				lastRec.setEnd(0);
				lastRec.partNumber = 0;
				lastRec.componentType = Constants::GAP;
				lastRec.gapLength = 0;
				lastRec.gapType = Constants::CONTIG;
				lastRec.isLinkage = false;
				RecList::const_iterator recIt;
				for (recIt = recs.begin(); recIt != recs.end(); ++recIt) {
					Record2& currRec = **recIt;
					if (currRec.isGap() and not currRec.isLinkage) {
						if (not lastRec.isGap()) {
							agpOutStream << lastRec;
							++recNum;
						}

						lastRec = currRec;
						lastRec.partNumber = recNum;
						agpOutStream << lastRec;			
						++recNum;
					} else if (lastRec.isGap()) {
						lastRec = currRec;
						lastRec.partNumber = recNum;
						lastRec.componentType = Constants::OTHER;
						lastRec.setSourceChrom(chrom + "." +
											   toString(chromRecNum));
						lastRec.setSourceStart(0);
						lastRec.setSourceEnd(lastRec.getLength());
						lastRec.setSourceStrand('+');
						++chromRecNum;
					} else {
						lastRec.setTargetEnd(currRec.getTargetEnd());
						lastRec.setSourceEnd(lastRec.getLength());
					}
				}

				if (not lastRec.getChrom().empty()) {
					agpOutStream << lastRec;
					++recNum;
				}
			}
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
