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

#include "polytope/Polytope.hh"
#include "polytope/formats/polymake/InputStream.hh"
#include "bio/alignment/AffineGapNWPairwiseConsensusAligner.hh"
#include "bio/alignment/AmbiguousDNAScoringMatrix.hh"
#include "bio/formats/fasta.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "filesystem.hh"
using namespace bio::alignment;
using namespace filesystem;
using namespace bio::formats;
using namespace polytope;
using util::string::toString;

#include "polyalign.hh"

typedef int T;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	std::string match = "0";
	std::string mismatch = "x";
	std::string space = "s";
	std::string gap = "g";
	std::string scoringMatrixFilename;
	std::string fastaFilename;
	std::string polytopeFilename;
	std::string outDirname = ".";

	// Set up option parser
	util::options::Parser parser("",
								 "Outputs a consensus alignment of the input "
								 "sequences for each normal cone of the input "
								 "polytope.");
	parser.addStoreOpt('m', "match", "match variable", match, "SCORE");
	parser.addStoreOpt('x', "mismatch", "mismatch variable", mismatch, "SCORE");
	parser.addStoreOpt('s', "space", "space variable", space, "SCORE");
	parser.addStoreOpt('g', "gap", "gap variable", gap, "SCORE");
	parser.addStoreOpt(0, "matrix", "scoring matrix filename",
					   scoringMatrixFilename, "FILENAME");

	parser.addStoreOpt('o', "out-dir",
					   "Directory into which to output vertex alignments",
					   outDirname, "DIR");
	parser.addStoreArg("fastaFile",
					   "FASTA file containing two sequences to be aligned",
					   fastaFilename);
	parser.addStoreArg("polytopeFile",
					   "POLYMAKE polytope file",
					   polytopeFilename);
	parser.parse(argv, argv + argc);

	try {
		// Read in sequences
		InputFileStream fastaFile(fastaFilename);
		fasta::InputStream fastaInputStream(fastaFile);
		fasta::Record seq1, seq2;
		fastaInputStream >> seq1 >> seq2;
		if (not fastaInputStream) {
			throw std::runtime_error("Could not read 2 sequences from "
									 "FASTA input");
		}

		// Read in polytope
		InputFileStream polytopeFile(polytopeFilename);
		formats::polymake::InputStream polymakeInputStream(polytopeFile);
		Polytope<T> p;
		polymakeInputStream >> p;

		// Set up scoring matrix, reading from file if specified
		VariableMatrix varMatrix;
		if (not scoringMatrixFilename.empty()) {
			InputFileStream scoringMatrixFile(scoringMatrixFilename);
			varMatrix.readMatrix(scoringMatrixFile);
		} else {
			varMatrix.setMatchScore(match);
			varMatrix.setMismatchScore(mismatch);
		}

		VariableSet varSet = makeVariableSet(varMatrix, space, gap);
		size_t numVars = getNumVars(varSet);
		
		// Check that dimensions of polytope correspond to number of params
		if (numVars != p.getAmbientDim()) {
			throw std::runtime_error("Invalid number of parameters " +
									 toString(numVars) + 
									 " for input polytope of dimension " +
									 toString(p.getAmbientDim()));
		}
		
		const Polytope<T>::FacetNormalList facetNormals = p.getFacetNormals();
		const Polytope<T>::FacetsThruVertexList facetsThruVertices = p.getFacetsThruVertices();

		Path outDir(outDirname);
		for (size_t i = 0; i < p.getNumVertices(); ++i) {
			// Construct a parameter ray that makes this vertex optimal
			Polytope<T>::Ray paramRay(p.getAmbientDim() + 1);
			typedef Polytope<T>::FacetsThruVertex::const_iterator FacetIndexIterator;
			for (FacetIndexIterator fi = facetsThruVertices[i].begin();
				 fi != facetsThruVertices[i].end(); ++fi) {
				paramRay += facetNormals[*fi];
			}

			// Initialize parameters according to the ray
			DNAScoringMatrix<int> scoreMatrix;
			int spaceScore, gapScore;
			setParameters(scoreMatrix, spaceScore, gapScore,
						  varMatrix, space, gap, varSet, paramRay);

			// Align sequences
			AffineGapNWPairwiseConsensusAligner<int>
				aligner(scoreMatrix, spaceScore, gapScore);
			PairwiseAlignment pa = aligner.align(seq1.sequence, seq2.sequence);

			// Output alignment
			fasta::Record r1(seq1.title, pa.seq1), r2(seq2.title, pa.seq2);
			OutputFileStream alignmentFile(outDir / (toString(i) + ".fa"));
			fasta::OutputStream fastaOutputStream(alignmentFile);
			fastaOutputStream << r1 << r2;
		}

	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
