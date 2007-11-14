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

#include "bio/agp/AGPBackwardMapper.hh"
#include "bio/genome/BasicInterval.hh"
#include "bio/gff/GFFRecord.hh"
#include "bio/gff/GFFInputStream.hh"
#include "util/string.hh"
#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;
using bio::genome::BasicInterval;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool remove = false;
	bool ignore = false;
	std::string segmentAttribute = "segment";
	std::string agpFilename;

	util::options::Parser parser("< gffInput", "");
	parser.addStoreTrueOpt('r', "remove",
						   "remove records that are not in a contig listed "
						   "in the AGP file",
						   remove);
	parser.addStoreTrueOpt('i', "ignore",
						   "ignore (do not transform) records that are not "
						   "in a contig listed in the AGP file",
						   ignore);
	parser.addStoreOpt('s', "segattr",
					   "name of attribute to add to mapped features that "
					   "broken up into multiple segments",
					   segmentAttribute, "ATTRIBUTE");
	parser.addStoreArg("agpFile", "", agpFilename);
	parser.parse(argv, argv + argc);

	try {
		InputFileStream agpFile(agpFilename);
		bio::agp::AGPBackwardMapper mapper(agpFile);
	
		// Construct GFF stream for fast reading
		std::cerr << "Untransforming GFF records...\n";
		bio::gff::GFFInputStream gffStream(std::cin);
		bio::gff::GFFRecord gffRec;
		while (gffStream >> gffRec) {
			std::vector<BasicInterval> mapped;
			mapper.map(gffRec.getInterval(), mapped);

			if (mapped.empty()) {
				if (remove) {
					continue;
				} else if (ignore) {
					std::cout << gffRec;
					continue;
				} else {
					std::cerr << "Warning: could not transform record:\n"
							  << gffRec << '\n';
					continue;
				}
			}

			for (size_t segNum = 0; segNum < mapped.size(); ++segNum) {
				bio::gff::GFFRecord untransformedRec = gffRec;
				untransformedRec.setInterval(mapped[segNum]);
				
				if (mapped.size() > 1) {
					untransformedRec.addAttribute(segmentAttribute).values.push_back(util::string::toString(segNum + 1));
				}

				std::cout << untransformedRec;
			}
		
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
