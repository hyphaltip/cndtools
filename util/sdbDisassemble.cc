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

struct AGPSourceStartSorter {
	bool operator()(const agp::Record* r1, const agp::Record* r2) const {
		return r1->getSourceStart() < r2->getSourceStart();
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
						   "use nib compression (DNA only) for disassembled DB",
						   compressed);
	parser.addStoreArg("sourceDB", "", sourceDBFilename);
	parser.addStoreArg("disassembledDB", "", outDBFilename);	
	parser.parse(argv, argv + argc);

	try {
		// Attempt to open source database
		bio::SDB::DB sourceDB;
		sourceDB.open(sourceDBFilename, true);
	
		// Attempt to open output database
		bio::SDB::DB outputDB;
		outputDB.open(outDBFilename, false, true);

		typedef std::vector<agp::Record*> RecList;
		typedef std::map<std::string, RecList> SourceRecsMap;
		SourceRecsMap chromRecs;
		
		// Read AGP records into map
		agp::InputStream agpStream(std::cin);
		agp::Record rec;
		while (agpStream >> rec) {
			if (not rec.isGap()) {
				agp::Record* sourceRec = new agp::Record(rec);
				chromRecs[rec.getSourceAccession()].push_back(sourceRec);
			}
		}
		
		// Construct each chromosome separately and insert into output DB
		SourceRecsMap::iterator it;
		for (it = chromRecs.begin(); it != chromRecs.end(); ++it) {
			const std::string& source = it->first;
			RecList& recs = it->second;
			
			// Order the AGP records by contig start coordinate
			std::sort(recs.begin(), recs.end(), AGPSourceStartSorter());
			
			// Determine the source length and make a new sequence
			std::string sourceSeq(recs.back()->getSourceEnd(), 'N');
			RecList::iterator recIt;
			bool sourceBuilt = false;
			for (recIt = recs.begin(); recIt != recs.end(); ++recIt) {
				std::string seq;
				try {
					seq = sourceDB.getSeq((*recIt)->getChrom(),
										  (*recIt)->getContigStart() - 1,
										  (*recIt)->getContigEnd(),
										  (*recIt)->getSourceOrientation());
				} catch (bio::SDB::NotFoundError& e) {
					std::cerr << "WARNING: " << e.what()
							  << ", when building record " << source << '\n';
					continue;
				}
				sourceBuilt = true;
				sourceSeq.replace((*recIt)->getSourceStart() - 1,
								 seq.size(),
								 seq);
			}
			if (sourceBuilt) {
				outputDB.putRec(source, sourceSeq, compressed);
			} else {
				std::cerr << "WARNING: Record " << source << " not built\n";
			}
		}
	} catch (std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		exit(EXIT_FAILURE);
	}				
		
	
	return EXIT_SUCCESS;
}
