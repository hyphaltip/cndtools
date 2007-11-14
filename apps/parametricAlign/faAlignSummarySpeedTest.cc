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
#include <string>
#include <vector>

#include "util/stl.hh"
#include "util/options.hh"
#include "filesystem.hh"
#include "bio/formats/fasta.hh"
#include "bio/alignment/DNAScoringMatrix.hh"
#include "boost/timer.hpp"

#include "AlignmentSummarizer.hh"

using bio::alphabet::DNA;
using namespace bio::alignment;
using namespace bio::formats::fasta;
using namespace filesystem;

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

	// Initialize options to defaults
	std::string match = "0";
	std::string mismatch = "x";
	std::string space = "s";
	std::string gap = "g";
	std::string scoringMatrixFilename;
	bool lexMin = false;
	bool linearMem = false;
	std::vector<int> paramValues;
	size_t iters = 100;

	// Parse command line
	util::options::Parser parser("< fastaInput", "");
	parser.addStoreTrueOpt(0, "lex-min", "Give lexicographically minimum summary", lexMin);
	parser.addStoreTrueOpt(0, "linear-mem", "Use linear space algorithm", linearMem);
	parser.addStoreOpt('m', "match", "match variable", match, "SYMBOL");
	parser.addStoreOpt('x', "mismatch", "mismatch variable", mismatch, "SYMBOL");
	parser.addStoreOpt('s', "space", "space variable", space, "SYMBOL");
	parser.addStoreOpt('g', "gap", "gap variable", gap, "SYMBOL");
	parser.addStoreOpt(0, "matrix", "symbolic scoring matrix filename",
					   scoringMatrixFilename, "FILENAME");
	parser.addStoreOpt(0, "iters", "iterations", iters);
	parser.addAppendArg("paramValue",
						"Value of a parameter in the model.  Parameter values "
						"must be given in alphabetical order of the parameter "
						"names",
						paramValues);
	parser.parse(argv, argv + argc);

	try {
		// Construct FASTA stream for fast reading
		InputStream fastaStream(std::cin);

		// Get first two records
		Record rec1, rec2;
		fastaStream >> rec1 >> rec2;
		if (not fastaStream) {
			throw std::runtime_error("Alignment requires two FASTA records");
		}

		// Check to make sure sequences only contain A, T, C, G
		checkSeq(rec1.sequence);
		checkSeq(rec2.sequence);
		
		DNAScoringMatrix<std::string> paramMatrix;
		if (not scoringMatrixFilename.empty()) {
			InputFileStream scoringMatrixFile(scoringMatrixFilename);
			paramMatrix.readMatrix(scoringMatrixFile);
		} else {
			paramMatrix.setMatchScore(match);
			paramMatrix.setMismatchScore(mismatch);
		}

		AlignmentSummarizer<int> summarizer(rec1.sequence, rec2.sequence,
											paramMatrix, space, gap,
											lexMin, linearMem);

		// Make sure the number of parameter values is the same as the
		// number of parameters
		paramValues.resize(summarizer.getNumParams());

		boost::timer timer;
		timer.restart();
		for (size_t i = 0; i < iters; ++i) {
			std::vector<size_t> s = summarizer.getSummary(paramValues);
		}
		std::cerr << timer.elapsed() << '\n';

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
