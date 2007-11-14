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

#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;

#include "BreakpointGraph.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options and arguments
	std::string breakpointFilename;
	
	// Set up option parser
	util::options::Parser parser("< roughHomologyMap > breakpointMultipleAlignments",
								 "Outputs segments of breakpoint regions that "
								 "should be multiply aligned according to the "
								 "given breakpoints");
	parser.addStoreArg("breakpointFile",
					   "File specifying breakpoint for each breakpoint segment",
					   breakpointFilename);
	parser.parse(argv, argv + argc);

	try {

		BreakpointGraph bg;
		
		bg.readHomologyMap(std::cin);
		InputFileStream breakpointFile(breakpointFilename);
		bg.setBreakpoints(breakpointFile);

		bg.writeAlignedSegments(std::cout);
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
