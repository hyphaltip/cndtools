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

#include "bio/alignment/AlphabetScoringMatrix.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "math/MaxPlus.hh"
#include "util/options.hh"
#include "filesystem.hh"
using bio::alignment::AlphabetScoringMatrix;

#include "AlignmentReader.hh"
#include "RecombinationScorer.hh"
#include "BestPathSemiRing.hh"

const bio::alphabet::Nucleotide GAPPED_DNA("ACGTNMRWSYKVHDB-",
										   "TGCANKYWSRMBDHV-");

void findBestPath(const MultipleAlignment& ma,
				  size_t recombinant,
				  int matchScore,
				  int mismatchScore,
				  int recombinationScore,
				  int norecombinationScore,
				  std::vector<size_t>& bestPath) {
	typedef BestPathSemiRing< math::MaxPlus<int> > SemiRing;
	typedef SemiRing::Element Element;
	math::MaxPlus<int> scoringSemiRing;
	SemiRing semiRing(scoringSemiRing);
	Element match(matchScore);
	Element mismatch(mismatchScore);
	std::vector<Element> initialScores(ma.getNumSeqs(),
									   semiRing.getMultiplicativeIdentity());
	AlphabetScoringMatrix<Element> emissionScores(GAPPED_DNA, match, mismatch);
	util::Matrix<Element> transitionScores(ma.getNumSeqs(),
										   ma.getNumSeqs());
	for (size_t i = 0; i < ma.getNumSeqs(); ++i) {
		for (size_t j = 0; j < ma.getNumSeqs(); ++j) {
			PathPtr p(new Path(j));
			transitionScores(i, j) = (i == j ? 
									  Element(norecombinationScore, p) :
									  Element(recombinationScore, p));
		}
	}

	RecombinationScorer<SemiRing>
		scorer(semiRing, initialScores, emissionScores, transitionScores);

	Element best = scorer.score(ma, recombinant);
	bestPath.clear();
	get_path_indices(best.path, bestPath);
}

const std::string USAGE = "";

const std::string DESCRIPTION =
"Given a multiple alignment, a selected sequence from the multiple "
"alignment to be a putative recombinant, and parameter settings, "
"returns *an* optimal path through the other sequences according to the "
"recombination model.";

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	size_t recombinant_num = 0;
	std::string recombinant_name = "";
	int match = 2;
	int mismatch = 0;
	int recombination = 0;
	int norecombination = 20;
	std::string alignment_filename;
	
	// Set up option parser
	util::options::Parser parser(USAGE, DESCRIPTION);
	parser.addStoreOpt(0, "recombinant",
					   "Number (0-based) of sequence to treat as recombinant",
					   recombinant_num, "INTEGER");
	parser.addStoreOpt(0, "recombinant-name",
					   "Name of sequence to treat as recombinant",
					   recombinant_name, "STRING");
	parser.addStoreOpt('m', "match", "Match score", match, "INTEGER");
	parser.addStoreOpt('x', "mismatch", "Mismatch score", mismatch, "INTEGER");
	parser.addStoreOpt('r', "recombination", "Recombination score",
					   recombination, "INTEGER");
	parser.addStoreOpt('n', "norecombination", "Norecombination score",
					   norecombination, "INTEGER");
	parser.addStoreArg("alignmentFile",
					   "Alignment file",
					   alignment_filename);
	parser.parse(argv, argv + argc);

	try {
		// Read in alignment
		filesystem::InputFileStream alignment_file(alignment_filename);
		AlignmentReader reader(alignment_file,
							   recombinant_name, recombinant_num);

		std::vector<size_t> bestPath;
		findBestPath(reader.getAlignment(), reader.getRecombinantNum(),
					 match, mismatch, recombination, norecombination,
					 bestPath);

		for (size_t i = 0; i < bestPath.size(); ++i) {
			if (i != 0) { std::cout << ' '; }
			std::cout << bestPath[i];
		}
		std::cout << '\n';

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
