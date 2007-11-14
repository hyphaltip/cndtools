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
#include "polytope/formats/polymake/OutputStream.hh"
#include "bio/formats/fasta.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "util/options.hh"
#include "filesystem.hh"
#include "util/string.hh"
using namespace filesystem;
using namespace bio::formats;
using namespace polytope;
using util::string::toString;

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

void outputPolytopeRow(PolytopeList& polytopes,
					   const std::string& prefix) {
	for (size_t i = 0; i < polytopes.size(); ++i) {
		OutputFileStream outFile(prefix + toString(i) + ".poly");
		formats::polymake::OutputStream polymakeOutputStream(outFile);
		polymakeOutputStream << polytopes[i];
	}
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	size_t numParams = 2;
	std::string qhullLogFile = "qhull.log";
	
	// Set up option parser
	util::options::Parser parser("< fastaInput > polymakeOutput",
								 DESCRIPTION);
	parser.addStoreOpt('p', "num-params",
					   "Number of parameters to use in alignment",
					   numParams, "NUM");
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

		// Align the sequences and get the polytope
		PolytopeList forwardMatch, forwardGap1, forwardGap2,
			backwardMatch, backwardGap1, backwardGap2;
		calcMiddleRowPolytopes(seq1.sequence,
							   seq2.sequence,
							   numParams,
							   forwardMatch,
							   forwardGap1,
							   forwardGap2,
							   backwardMatch,
							   backwardGap1,
							   backwardGap2);
		outputPolytopeRow(forwardMatch, "fh");
		outputPolytopeRow(forwardGap1, "fd");
		outputPolytopeRow(forwardGap2, "fi");
		outputPolytopeRow(backwardMatch, "bh");
		outputPolytopeRow(backwardGap1, "bd");
		outputPolytopeRow(backwardGap2, "bi");

	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
