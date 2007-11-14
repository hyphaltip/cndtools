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
#include "bio/agp/AGPBackwardMapper.hh"
#include "util/io/tabdelim.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "util/stl.hh"
#include "filesystem.hh"
using namespace filesystem;
using util::string::toString;
using namespace bio;
using bio::genome::BasicInterval;

std::vector<size_t> parse_string_indices(const std::string& s) {
	std::vector<std::string> string_indices;
	util::string::split(s, std::back_inserter(string_indices), ",");
	std::vector<size_t> indices;
	util::string::Converter<size_t> toNum;
	for (size_t i = 0; i < string_indices.size(); ++i) {
		indices.push_back(toNum(string_indices[i]) - 1);
	}
	return indices;
}

BasicInterval make_interval(const std::vector<std::string>& fields,
							const std::vector<size_t> indices) {
	util::string::Converter<genome::Position> toPosition;
	return BasicInterval(fields[indices[0]],
						 toPosition(fields[indices[1]]),
						 toPosition(fields[indices[2]]),
						 fields[indices[3]].at(0));
}

void transform_interval(std::vector<std::string>& fields,
						const std::vector<size_t> indices,
						const BasicInterval& interval) {
	fields[indices[0]] = interval.getChrom();
	fields[indices[1]] = toString(interval.getStart());
	fields[indices[2]] = toString(interval.getEnd());
	fields[indices[3]] = toString(interval.getStrand());
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::string stringIndices = "1,2,3,4";
	std::string agpFilename;
	bool backward = false;

	util::options::Parser parser("< tabdelimInput", "");
	parser.addStoreOpt('f', "fields",
					   "Fields containing the genome interval to transform",
					   stringIndices);
	parser.addStoreTrueOpt('b', "backward",
						   "Transform backwards (chromosomes to contigs)",
						   backward);
	parser.addStoreArg("agpFile", "", agpFilename);
	parser.parse(argv, argv + argc);

	try {
		std::vector<size_t> indices(parse_string_indices(stringIndices));
		if (indices.size() != 4) {
			throw std::runtime_error("Must specify 4 fields for interval");
		}
		
		InputFileStream agpFile(agpFilename);

		bio::genome::IntervalMapper* mapper;
		if (backward) {
			mapper = new bio::agp::AGPBackwardMapper(agpFile);
		} else {
			mapper = new bio::agp::AGPForwardMapper(agpFile);
		}

		util::io::tabdelim::InputStream fieldsInputStream(std::cin);
		util::io::tabdelim::OutputStream fieldsOutputStream(std::cout);
		std::vector<std::string> fields;
		std::vector<BasicInterval> mapping;

		while (fieldsInputStream >> fields) {
			mapping.clear();
			BasicInterval interval = make_interval(fields, indices);
			mapper->map(interval, mapping);
			if (mapping.size() != 1) {
				util::stl::print_elements(std::cerr, mapping);
				std::cerr << mapping.size() << '\n';
				throw std::runtime_error("Could not transform interval: " +
										 toString(interval));
			}
			transform_interval(fields, indices, mapping.front());
			fieldsOutputStream << fields;
		}

		delete mapper;

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
