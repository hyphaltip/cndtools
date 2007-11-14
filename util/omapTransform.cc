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
#include <fstream>

#include "bio/formats/agp/Record.hh"
#include "bio/agp/AGPForwardMapper.hh"
#include "bio/homologymap/Map.hh"
#include "util/options.hh"
#include "util/string.hh"
using util::string::toString;

using bio::genome::BasicInterval;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::vector<std::string> agpFilenames;

	util::options::Parser parser("< orthologyMap", "");
	parser.addAppendArg("agpFilename", "", agpFilenames);
	parser.parse(argv, argv + argc);

	try {
		std::vector<bio::agp::AGPForwardMapper*> mappers;
		for (size_t i = 0; i < agpFilenames.size(); ++i) {
			std::ifstream agpFile(agpFilenames[i].c_str());
			mappers.push_back(new bio::agp::AGPForwardMapper(agpFile));
		}

		bio::homologymap::Segment seg;
		while (std::cin >> seg) {
			assert(seg.intervals.size() == mappers.size());
			for (size_t i = 0; i < seg.intervals.size(); ++i) {
				if (not seg.hasGenome(i)) {
					continue;
				}
				std::vector<BasicInterval> mappings;
				mappers[i]->map(*seg.intervals[i], mappings);
				if (mappings.size() == 0) {
					throw std::runtime_error("Mapping for interval: " +
											 toString(*seg.intervals[i]) +
											 " does not exist");
				} else if (mappings.size() > 1) {
					throw std::runtime_error("Mapping for interval: " +
											 toString(*seg.intervals[i]) +
											 " is larger than 1");
				}					
				delete seg.intervals[i];
				seg.intervals[i] = new BasicInterval(mappings[0]);
			}
			std::cout << seg;
			for (size_t i = 0; i < seg.intervals.size(); ++i) {
				if (not seg.hasGenome(i)) {
					continue;
				}
				delete seg.intervals[i];
			}
		}
	} catch (std::runtime_error& e) {
		std::cerr << "ERROR: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
