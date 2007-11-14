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

#include "bio/agp/AGPForwardMapper.hh"
#include "bio/gff/GFFRecord.hh"
#include "bio/gff/GFFInputStream.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "filesystem.hh"
using util::string::toString;
using bio::genome::BasicInterval;
using namespace filesystem;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool remove = false;
	bool ignore = false;
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
	parser.addStoreArg("agpFile", "", agpFilename);
	parser.parse(argv, argv + argc);

	try {
		InputFileStream agpFile(agpFilename);
		bio::agp::AGPForwardMapper mapper(agpFile);
		
		// Construct GFF stream for fast reading
		bio::gff::GFFInputStream gffStream(std::cin);
		bio::gff::GFFRecord gffRec;
		std::vector<BasicInterval> mapping;	
		while (gffStream >> gffRec) {
			mapping.clear();
			mapper.map(gffRec.getInterval(), mapping);
			if (mapping.empty()) {
				if (remove) {
					continue;
				} else if (ignore) {
					std::cout << gffRec;
					continue;
				} else {
					throw std::runtime_error(gffRec.getSeqname() +
											 " not listed in AGP file");
				}
			} else if (mapping.size() > 1) {
				std::cerr << "GFF record spans multiple intervals:\n"
						  << gffRec
						  << mapping[0] << "\n"
						  << mapping[1] << '\n';
				continue;
			}
			
			gffRec.setSeqname(mapping[0].getChrom());
			gffRec.setStart(mapping[0].getStart() + 1);
			gffRec.setEnd(mapping[0].getEnd());
			if (gffRec.hasStrand()) {
				gffRec.setStrand(mapping[0].getStrand());
			}
			std::cout << gffRec;
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
