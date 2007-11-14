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

#include "polytope/Polytope.hh"
#include "polytope/formats/polymake/InputStream.hh"
#include "polytope/formats/polymake/OutputStream.hh"
#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;
using namespace polytope;

typedef int T;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	std::vector<std::string> polytopeFilenames;
	
	// Set up option parser
	util::options::Parser parser("> polytopeOutput",
								 "Outputs the Minkowski sum of the input "
								 "polytopes");
	parser.addAppendArg("polytopeFile",
						"POLYMAKE polytope file",
						polytopeFilenames);
	parser.parse(argv, argv + argc);

	try {
		Polytope<T> product;
		bool first = true;
		
		for (size_t i = 0; i < polytopeFilenames.size(); ++i) {
			InputFileStream polytopeFile(polytopeFilenames[i]);
			formats::polymake::InputStream polymakeInputStream(polytopeFile);
			Polytope<T> p;
			polymakeInputStream >> p;
			if (first) {
				product = p;
				first = false;
			} else {
				product *= p;
			}
		}

		formats::polymake::OutputStream polymakeOutputStream(std::cout);
		polymakeOutputStream << product;
		
	}  catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
