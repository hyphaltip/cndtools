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
#include <vector>
#include <stdexcept>

#include "util/matrix.hh"
#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;

#include "BreakpointGraph.hh"
#include "DistanceCalculator.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options and arguments
	std::string genomeDir = ".";
	std::string homologyMapFilename;
	std::string treeFilename;
	std::string segmentsFilename = "segments";
	std::string edgesFilename = "edges";
	bool removeColinearEdges = false;
	
	// Set up option parser
	util::options::Parser parser("",
								 "Finds breakpoints in between homologous "
								 "segments in multiple genomes");
	parser.addStoreOpt(0,
					   "genome-dir",
					   "directory containing genome sequences in sdb format",
					   genomeDir,
					   "DIR");
	parser.addStoreOpt(0,
					   "segments",
					   "file to output breakpoint segment coordinates",
					   segmentsFilename,
					   "FILE");
	parser.addStoreTrueOpt(0,
						   "remove-colinear",
						   "Remove edges in the graph that connect colinear "
						   "segments",
						   removeColinearEdges);
	parser.addStoreOpt(0,
					   "edges",
					   "file to output edges of the breakpoint graph",
					   edgesFilename,
					   "FILE");
	parser.addStoreArg("homologyMapFile",
					   "Homology map specifying the homlogous segments "
					   "between which breakpoints will be found",
					   homologyMapFilename);
	parser.addStoreArg("treeFile",
					   "Newick-formatted tree containing the species in the "
					   "homology map",
					   treeFilename);
	parser.parse(argv, argv + argc);

	try {
		Genome::setDBDir(genomeDir);

		InputFileStream treeFile(treeFilename);
		DistanceCalculator dc(treeFile);
		util::Matrix<double> distances = dc.getDistances();

		InputFileStream homologyMapFile(homologyMapFilename);

		BreakpointGraph bg;
		bg.readHomologyMap(homologyMapFile);
		bg.addAllEdges(distances, removeColinearEdges);

		OutputFileStream segmentsFile(segmentsFilename);
		bg.writeBreakpointSegments(segmentsFile);

		OutputFileStream edgesFile(edgesFilename);
		bg.writeMSTEdges(edgesFile);
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
