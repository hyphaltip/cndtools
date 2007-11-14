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
#include <map>
#include <vector>

#include "bio/formats/agp/Record.hh"
#include "bio/formats/agp/InputStream.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "bio/sdb.hh"
#include "util/options.hh"
using namespace bio::formats;

struct AGPContigStartSorter {
	bool operator()(const agp::Record* r1, const agp::Record* r2) const {
		return r1->getContigStart() < r2->getContigStart();
	}
};

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool compressed = false;
	std::string sourceDBFilename;
	std::string outDBFilename;
	
	util::options::Parser parser("< agpInput", "");
	parser.addStoreTrueOpt('c', "compressed",
						   "use nib compression (DNA only) for assembled DB",
						   compressed);
	parser.addStoreArg("sourceDB", "", sourceDBFilename);
	parser.addStoreArg("assembledDB", "", outDBFilename);	
	parser.parse(argv, argv + argc);

	try {
		// Attempt to open source database
		bio::SDB::DB sourceDB;
		sourceDB.open(sourceDBFilename, true);
	
		// Attempt to open output database
		bio::SDB::DB outputDB;
		outputDB.open(outDBFilename, false, true);

		typedef std::vector<agp::Record*> RecList;
		typedef std::map<std::string, RecList> ChromRecsMap;
		ChromRecsMap chromRecs;
		
		// Read AGP records into map
		agp::InputStream agpStream(std::cin);
		agp::Record rec;
		while (agpStream >> rec) {
			chromRecs[rec.getChrom()].push_back(new agp::Record(rec));
		}
		
		// Construct each chromosome separately and insert into output DB
		ChromRecsMap::iterator it;
		for (it = chromRecs.begin(); it != chromRecs.end(); ++it) {
			const std::string& chrom = it->first;
			RecList& recs = it->second;
			
			// Order the AGP records by contig start coordinate
			std::sort(recs.begin(), recs.end(), AGPContigStartSorter());
			
			// Determine the chrom length and make a new sequence
			std::string chromSeq(recs.back()->getContigEnd(), 'N');
			RecList::iterator recIt;
			for (recIt = recs.begin(); recIt != recs.end(); ++recIt) {
				if ((*recIt)->isGap()) {
					continue;
				}

				std::string seq;
				seq = sourceDB.getSeq((*recIt)->getSourceAccession(),
									  (*recIt)->getSourceStart() - 1,
									  (*recIt)->getSourceEnd(),
									  (*recIt)->getSourceOrientation());

				chromSeq.replace((*recIt)->getContigStart() - 1, seq.size(), seq);
			}
			outputDB.putRec(chrom, chromSeq, compressed);
		}
	} catch (std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		exit(EXIT_FAILURE);
	}				
		
	
	return EXIT_SUCCESS;
}
