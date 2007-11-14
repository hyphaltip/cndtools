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
#include "filesystem.hh"
#include "bio/formats/fasta.hh"
#include "bio/alignment/AffineGapNWPairwiseConsensusAligner.hh"
#include "bio/alignment/AffineGapNWPairwiseAligner.hh"
#include "bio/alignment/AffineGapNWFastPairwiseAligner.hh"
#include "bio/alignment/AffineGapNWFastPairwiseAligner2.hh"
#include "bio/alignment/AffineGapNWFastPairwiseAligner3.hh"
#include "bio/alignment/AmbiguousDNAScoringMatrix.hh"

using namespace bio::alignment;
using namespace bio::formats::fasta;
using namespace filesystem;

void outputAlignment(OutputStream& strm,
					 const PairwiseAlignment& pa,
					 const std::string& title1,
					 const std::string& title2) {
	Record a1(title1, pa.seq1), a2(title2, pa.seq2);
	strm << a1 << a2;
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	int match = 100;
	int mismatch = -100;
	int space = -30;
	int gap = -400;
	std::string scoringMatrixFilename;
	bool all = false;
	bool consensus = false;
	bool fast = false;
	bool fast2 = false;
	bool fast3 = false;

	// Parse command line
	util::options::Parser parser("< fastaInput > fastaOutput", "");
	parser.addStoreOpt('m', "match", "match score", match, "SCORE");
	parser.addStoreOpt('x', "mismatch", "mismatch score", mismatch, "SCORE");
	parser.addStoreOpt('s', "space", "space score", space, "SCORE");
	parser.addStoreOpt('g', "gap", "gap score", gap, "SCORE");
	parser.addStoreOpt(0, "matrix", "scoring matrix filename",
					   scoringMatrixFilename, "FILENAME");
	parser.addStoreTrueOpt(0, "fast", "Use fast (quadratic space) algorithm",
						   fast);
	parser.addStoreTrueOpt(0, "fast2", "use fast2", fast2);
	parser.addStoreTrueOpt(0, "fast3", "use fast3", fast3);
	parser.addStoreTrueOpt('a', "all", "give all optimal alignments", all);
	parser.addStoreTrueOpt('c', "consensus", "give consensus alignment",
						   consensus);
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

		DNAScoringMatrix<int> dnaMatrix;
		if (not scoringMatrixFilename.empty()) {
			InputFileStream scoringMatrixFile(scoringMatrixFilename);
			dnaMatrix.readMatrix(scoringMatrixFile);
		} else {
			dnaMatrix.setMatchScore(match);
			dnaMatrix.setMismatchScore(mismatch);
		}
		AmbiguousDNAScoringMatrix<int> matrix(dnaMatrix);

		OutputStream fastaOutStream(std::cout);
		
		if (all) {
			AffineGapNWPairwiseAligner<int> aligner(matrix, space, gap);

			typedef std::vector<PairwiseAlignment> AlignmentList;
			AlignmentList alignments;
			aligner.align(rec1.sequence, rec2.sequence, alignments);
						
			for (AlignmentList::const_iterator it = alignments.begin();
				 it != alignments.end(); ++it) {
				outputAlignment(fastaOutStream, *it, rec1.title, rec2.title);
			}
		} else if (consensus) {
			AffineGapNWPairwiseConsensusAligner<int> aligner(matrix, space, gap);
			PairwiseAlignment pa = aligner.align(rec1.sequence, rec2.sequence);
			outputAlignment(fastaOutStream, pa, rec1.title, rec2.title);
		} else {
			PairwiseAligner* aligner;

			if (fast) {
				aligner = new AffineGapNWFastPairwiseAligner<int>(matrix, space, gap);
			} else if (fast2) {
				aligner = new AffineGapNWFastPairwiseAligner2<int>(matrix, space, gap);
			} else if (fast3) {
				aligner = new AffineGapNWFastPairwiseAligner3<int>(matrix, space, gap);
			} else {
				aligner = new AffineGapNWPairwiseAligner<int>(matrix, space, gap);
			}

			PairwiseAlignment pa = aligner->align(rec1.sequence, rec2.sequence);
			outputAlignment(fastaOutStream, pa, rec1.title, rec2.title);
			delete aligner;
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
