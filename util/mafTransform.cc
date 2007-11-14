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

#include "bio/formats/maf.hh"
#include "bio/agp/AGPForwardMapper.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "filesystem.hh"
using namespace filesystem;
using util::string::toString;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::string agpFilename;
	
	// Parse command line
	util::options::Parser parser("< mafInput", "Transform MAF input");
	parser.addStoreArg("agpFile", "", agpFilename);
	parser.parse(argv, argv + argc);

	try {
		InputFileStream agpFile(agpFilename);
		bio::agp::AGPForwardMapper mapper(agpFile);

		bio::formats::maf::InputStream inputStream(std::cin);
		bio::formats::maf::OutputStream outputStream(std::cout,
													 inputStream.getHeader());
		
		std::vector<bio::genome::BasicInterval> mapping;

		bio::formats::maf::Record rec;
		while (inputStream >> rec) {
			for (size_t i = 0; i < rec.sequences.size(); ++i) {
				bio::formats::maf::Sequence& seq = rec.sequences[i];
				mapping.clear();
				mapper.map(seq, mapping);
				if (mapping.size() != 1) {
					throw std::runtime_error("Could not map interval: " +
											 toString(seq));
				}
				seq.srcSize = mapper.getChromSize(mapping.front().getChrom());
				seq.setInterval(mapping.front());
			}

			outputStream << rec;
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
