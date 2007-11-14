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
#include <map>
#include <set>

#include "polytope/Polytope.hh"
#include "polytope/formats/polymake/OutputStream.hh"
#include "bio/formats/fasta.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "bio/alignment/DNAScoringMatrix.hh"
#include "util/options.hh"
#include "filesystem.hh"
using namespace filesystem;
using namespace bio::formats;
using namespace polytope;
using bio::alphabet::DNA;
using bio::alignment::DNAScoringMatrix;

#include "polyalign.hh"

const std::string DESCRIPTION = 
"Calculates the newton polytope for the polynomial representing all possible "
"alignments of the two sequences given as arguments.";


// Check to make sure the sequence in REC contains only DNA characters
void checkSeq(const bio::formats::fasta::Record& rec) {
	if (!bio::alphabet::DNA.isOn(rec.sequence)) {
		throw std::runtime_error("Sequence has characters other than "
								 "A, T, C, G: " + rec.title);
	}
}				   

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	std::string match = "0";
	std::string mismatch = "x";
	std::string space = "s";
	std::string gap = "g";
	std::string scoringMatrixFilename;
	std::string qhullLogFile = "qhull.log";
	
	// Set up option parser
	util::options::Parser parser("< fastaInput > polymakeOutput",
								 DESCRIPTION);
	parser.addStoreOpt('m', "match", "match variable", match, "SCORE");
	parser.addStoreOpt('x', "mismatch", "mismatch variable", mismatch, "SCORE");
	parser.addStoreOpt('s', "space", "space variable", space, "SCORE");
	parser.addStoreOpt('g', "gap", "gap variable", gap, "SCORE");
	parser.addStoreOpt(0, "matrix", "scoring matrix filename",
					   scoringMatrixFilename, "FILENAME");
	parser.addStoreOpt(0, "qhull-log",
					   "File for qhull warning/error messages",
					   qhullLogFile, "FILE");
	parser.parse(argv, argv + argc);

	try {
		Polytope<int>::setQhullLogFile(qhullLogFile);
		
		fasta::Record seq1, seq2;
		fasta::InputStream fastaInputStream(std::cin);
		fastaInputStream >> seq1 >> seq2;
		if (not fastaInputStream) {
			throw std::runtime_error("Could not read 2 sequences from "
									 "FASTA input");
		}
				
		// Check to make sure sequences only contain A, T, C, G
		checkSeq(seq1);
		checkSeq(seq2);

		VariableMatrix varMatrix;
		if (not scoringMatrixFilename.empty()) {
			InputFileStream scoringMatrixFile(scoringMatrixFilename);
			varMatrix.readMatrix(scoringMatrixFile);
		} else {
			varMatrix.setMatchScore(match);
			varMatrix.setMismatchScore(mismatch);
		}

		VariableSet varSet = makeVariableSet(varMatrix, space, gap);
		VariablePolytopeMap varMap = makeVariablePolytopeMap(varSet);
		PolytopeMatrix matrix = makePolytopeMatrix(varMatrix, varMap);
		size_t numVars = getNumVars(varSet);
		
		// Align the sequences and get the polytope
		IntegerPolytope p = calculatePolytope(seq1.sequence,
											  seq2.sequence,
											  matrix,
											  varMap[space],
											  varMap[gap],
											  numVars);

		formats::polymake::OutputStream polymakeOutputStream(std::cout);

		polymakeOutputStream << p;

	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
