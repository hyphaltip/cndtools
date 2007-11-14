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

#include "util/options.hh"
#include "util/string.hh"
#include "filesystem.hh"
using namespace filesystem;
using util::string::toString;

#include "AlignmentReader.hh"
#include "BestStateSetLister.hh"
#include "OrderedPolygon.hh"

const std::string USAGE = "";

const std::string DESCRIPTION =
"Given a multiple alignment, a selected sequence from the multiple "
"alignment to be a putative recombinant, and the recombination netwon "
"polygon for these inputs, outputs the normal cone for each vertex of the "
"poltyope as well as the optimal sequence subsets for that cone.";

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	size_t recombinant_num = 0;
	std::string recombinant_name = "";
	std::string alignment_filename;
	std::string polygon_filename;

	// Set up option parser
	util::options::Parser parser(USAGE, DESCRIPTION);
	parser.addStoreOpt(0, "recombinant",
					   "Number (0-based) of sequence to treat as recombinant",
					   recombinant_num, "INTEGER");
	parser.addStoreOpt(0, "recombinant-name",
					   "Name of sequence to treat as recombinant",
					   recombinant_name, "STRING");
	parser.addStoreArg("alignmentFile",
					   "Alignment file",
					   alignment_filename);
	parser.addStoreArg("polygonFile",
					   "POLYMAKE polygon file",
					   polygon_filename);
	parser.parse(argv, argv + argc);

	try {
		// Read in alignment
		InputFileStream alignment_file(alignment_filename);
		AlignmentReader reader(alignment_file,
							   recombinant_name, recombinant_num);
		
		// Read in polygon
		InputFileStream polygon_file(polygon_filename);
		OrderedPolygon<int> p;
		polygon_file >> p;

		for (size_t i = 0; i < p.getNumVertices(); ++i) {
			Vector<int> normal = p.vertexNormal(i);
			BestStateSetLister lister(normal[0],  0, 0, normal[1]);
			StateSetList bestStateSetList;
			lister.getBestStateSetList(reader.getAlignment(),
									   reader.getRecombinantNum(),
									   bestStateSetList);

			for (size_t i = 0; i < bestStateSetList.size(); ++i) {
				StateSet ss = bestStateSetList[i];
				size_t first = ss.find_first();
				for (size_t state = first; state < ss.size(); ++state) {
					if (not ss[state]) { continue; }
					if (state != first) { std::cout << ','; }
					std::cout << state;
				}
				std::cout << '\t';
			}

			std::cout << '\n';
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
