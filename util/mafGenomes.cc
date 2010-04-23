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
#include "boost/unordered_set.hpp"
#include "util/options.hh"
using namespace bio::formats;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	
	// Parse command line
	util::options::Parser parser("< mafInput", "List names of genomes in MAF");
	parser.parse(argv, argv + argc);

	try {
		maf::InputStream input_stream(std::cin);
        boost::unordered_set<std::string> genomes;
		
		maf::Record rec;
		while (input_stream >> rec) {
			for (size_t i = 0; i < rec.sequences.size(); ++i) {
				genomes.insert(rec.sequences[i].getGenomeAndChrom().first);
			}
		}

		util::stl::print_elements(std::cout, genomes, "\n");
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
