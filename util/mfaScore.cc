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

#include "bio/alignment/BasicNamedMultipleAlignment.hh"
#include "bio/alignment/AffineGapPairwiseAlignmentScorer.hh"
#include "bio/alignment/SumOfPairsScorer.hh"
#include "bio/alignment/AlphabetScoringMatrix.hh"
#include "bio/alignment/AmbiguousDNAScoringMatrix.hh"
#include "bio/alphabet/Alphabet.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "bio/alphabet/AminoAcid.hh"
#include "bio/formats/fasta.hh"
#include "filesystem.hh"
#include "util/options.hh"

using namespace bio::alignment;
using namespace bio::formats;
using namespace filesystem;

typedef long long Score;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);
	
	// Initialize options to defaults
	bool showFilenames = false;
	bool showTotal = false;
	bool onlyTotal = false;
	Score match = 100;
	Score mismatch = -100;
	Score space = -30;
	Score gap = -400;
	bool protein = false;
	std::string scoringMatrixFilename;
	std::vector<std::string> filenames;

	// Parse command line
	util::options::Parser parser("", "Score multiple alignments");
	parser.addStoreOpt('m', "match", "match score", match, "SCORE");
	parser.addStoreOpt('x', "mismatch", "mismatch score", mismatch, "SCORE");
	parser.addStoreOpt('s', "space", "space score", space, "SCORE");
	parser.addStoreOpt('g', "gap", "gap score", gap, "SCORE");
	parser.addStoreOpt(0, "matrix", "scoring matrix filename",
					   scoringMatrixFilename, "FILENAME");
	parser.addStoreTrueOpt(0, "prot", "score protein alignment", protein);
	parser.addStoreFalseOpt(0, "dna", "score dna alignment (default)", protein);
	parser.addStoreTrueOpt(0, "show-filenames",
						   "show the filename for each multiple alignment score",
						   showFilenames);
	parser.addStoreTrueOpt(0, "show-total",
						   "show the total score as the last line", showTotal);
	parser.addStoreTrueOpt(0, "only-total",
						   "only show the total score", onlyTotal);
	parser.addAppendArg("mfaFile",
						"multiple alignment file",
						filenames);
	
	parser.parse(argv, argv + argc);

	try {
		ScoringMatrix<Score>* matrix;

		if (protein) {
			AlphabetScoringMatrix<Score>* aaMatrix = new AlphabetScoringMatrix<Score>(bio::alphabet::UnambiguousAminoAcid);
			if (not scoringMatrixFilename.empty()) {
				std::ifstream scoringMatrixFile(scoringMatrixFilename.c_str());
				aaMatrix->readMatrix(scoringMatrixFile);
			} else {
				aaMatrix->setMatchScore(match);
				aaMatrix->setMismatchScore(mismatch);
			}
			matrix = aaMatrix;
		} else {
			DNAScoringMatrix<Score> dnaMatrix;
			if (not scoringMatrixFilename.empty()) {
				std::ifstream scoringMatrixFile(scoringMatrixFilename.c_str());
				dnaMatrix.readMatrix(scoringMatrixFile);
			} else {
				dnaMatrix.setMatchScore(match);
				dnaMatrix.setMismatchScore(mismatch);
			}
			matrix = new AmbiguousDNAScoringMatrix<Score>(dnaMatrix);
		}

		AffineGapPairwiseAlignmentScorer<Score> pairwiseScorer(*matrix,
															   space, gap);
		SumOfPairsScorer<Score> scorer(pairwiseScorer);

		Score totalScore = 0;
		for (size_t i = 0; i < filenames.size(); ++i) {
			InputFileStream file(filenames[i]);
			fasta::InputStream fastaStream(file);
			BasicNamedMultipleAlignment ma;
			fastaStream >> ma;
			Score score = scorer.score(ma);
			totalScore += score;
			if (not onlyTotal) {
				if (showFilenames) {
					std::cout << filenames[i] << '\t';
				}
				std::cout << score << std::endl;
			}
		}

		if (showTotal or onlyTotal) {
			if (showFilenames) {
				std::cout << "Total" << '\t';
			}
			std::cout << totalScore << std::endl;
		}
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
