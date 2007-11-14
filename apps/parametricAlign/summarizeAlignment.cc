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
#include "bio/alignment/PairwiseAlignment.hh"

#include "AlignmentSummarizer.hh"

using bio::alphabet::DNA;
using namespace bio::alignment;
using namespace bio::formats::fasta;
using namespace filesystem;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::string match = "0";
	std::string mismatch = "x";
	std::string space = "s";
	std::string gap = "g";
	std::string scoringMatrixFilename;

	// Parse command line
	util::options::Parser parser("< alignmentInput", "");
	parser.addStoreOpt('m', "match", "match variable", match, "SYMBOL");
	parser.addStoreOpt('x', "mismatch", "mismatch variable", mismatch, "SYMBOL");
	parser.addStoreOpt('s', "space", "space variable", space, "SYMBOL");
	parser.addStoreOpt('g', "gap", "gap variable", gap, "SYMBOL");
	parser.addStoreOpt(0, "matrix", "symbolic scoring matrix filename",
					   scoringMatrixFilename, "FILENAME");
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

		DNAScoringMatrix<std::string> paramMatrix;
		if (not scoringMatrixFilename.empty()) {
			InputFileStream scoringMatrixFile(scoringMatrixFilename);
			paramMatrix.readMatrix(scoringMatrixFile);
		} else {
			paramMatrix.setMatchScore(match);
			paramMatrix.setMismatchScore(mismatch);
		}

		AlignmentSummarizer<int> summarizer("", "", paramMatrix, space, gap);

		PairwiseAlignment pa(rec1.sequence, rec2.sequence);
		
		std::vector<size_t> s = summarizer.getSummary(pa);

		util::stl::print_elements(std::cout, s, "\t");

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
