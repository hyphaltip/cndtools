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
#include <stdexcept>

#include "bio/alignment/DNAScoringMatrix.hh"
#include "bio/alignment/AmbiguousDNAScoringMatrix.hh"
#include "util/matrix.hh"
#include "util/options.hh"
#include "filesystem.hh"
using bio::alignment::DNAScoringMatrix;
using bio::alignment::AmbiguousDNAScoringMatrix;
using namespace filesystem;

#include "BreakpointGraph.hh"
#include "DistanceCalculator.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options and arguments
	std::string genomeDir = ".";
	size_t maxSegmentPositions = 100;
	std::string homologyMapFilename;
	std::string treeFilename;
	std::string edgesFilename;
	std::string alignmentsDir;
	std::string alignmentFilename = "mavid.mfa";
	Score match = 100;
	Score mismatch = -100;
	Score space = -30;
	Score gap = -400;
	std::string scoringMatrixFilename;
	
	// Set up option parser
	util::options::Parser parser("",
								 "Finds breakpoints in between homologous "
								 "segments in multiple genomes");
	parser.addStoreOpt(0,
					   "genome-dir",
					   "directory containing genome sequences in sdb format",
					   genomeDir,
					   "DIR");
	parser.addStoreOpt('m', "match", "match score", match, "SCORE");
	parser.addStoreOpt('x', "mismatch", "mismatch score", mismatch, "SCORE");
	parser.addStoreOpt('s', "space", "space score", space, "SCORE");
	parser.addStoreOpt('g', "gap", "gap score", gap, "SCORE");
	parser.addStoreOpt(0, "matrix", "scoring matrix filename",
					   scoringMatrixFilename, "FILENAME");
	parser.addStoreOpt('r',
					   "resolution",
					   "number of breakpoint positions to consider within "
					   "each breakpoint segment on each iteration",
					   maxSegmentPositions,
					   "NUM");
	parser.addStoreOpt(0, "alignmentFile",
					   "Name of alignment file in each edge directory",
					   alignmentFilename, "STRING");
	parser.addStoreArg("homologyMapFile",
					   "Homology map specifying the homlogous segments "
					   "between which breakpoints will be found",
					   homologyMapFilename);
	parser.addStoreArg("treeFile",
					   "Newick-formatted tree containing the species in the "
					   "homology map",
					   treeFilename);
	parser.addStoreArg("edgesFile",
					   "File containing the edges of the breakpoint graph",
					   edgesFilename);
	parser.addStoreArg("alignmentsDir",
					   "Directory containing alignments for the edges",
					   alignmentsDir);
	parser.parse(argv, argv + argc);

	try {
		Genome::setDBDir(genomeDir);
		AlignedSegments::setAlignmentsDir(alignmentsDir);
		AlignedSegments::setAlignmentFilename(alignmentFilename);

		DNAScoringMatrix<Score> dnaMatrix;
		if (not scoringMatrixFilename.empty()) {
			InputFileStream scoringMatrixFile(scoringMatrixFilename);
			dnaMatrix.readMatrix(scoringMatrixFile);
		} else {
			dnaMatrix.setMatchScore(match);
			dnaMatrix.setMismatchScore(mismatch);
		}
		AmbiguousDNAScoringMatrix<Score> matrix(dnaMatrix);
		
		AlignedSegments::setScores(matrix, space, gap);

		InputFileStream treeFile(treeFilename);
		DistanceCalculator dc(treeFile);
		util::Matrix<double> distances = dc.getDistances();

		InputFileStream homologyMapFile(homologyMapFilename);

		BreakpointGraph bg;
		
		bg.readHomologyMap(homologyMapFile);

		InputFileStream edgesFile(edgesFilename);
		bg.addSomeEdges(edgesFile, distances);

		bg.calcMidpointScore(distances);
		
		bg.calcBestBreakpoints(distances, maxSegmentPositions);
		bg.writeBreakpoints(std::cout);
		
	} catch (const std::runtime_error& e) {
		std::cerr << "\nError: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
