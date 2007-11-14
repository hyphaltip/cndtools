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

#include "boost/random.hpp"
#include "boost/timer.hpp"

#include "util/options.hh"
#include "filesystem.hh"

#include "polyalign.hh"

using namespace filesystem;

const std::string DNA = "ACGT";

template<typename RNG>
std::string randomDNA(size_t length, RNG& rng) {
	std::string seq(length, ' ');
	for (size_t i = 0; i < length; ++i) {
		seq[i] = DNA[rng()];
	}
	return seq;
}

template<typename RNG>
std::string simpleRepeatDNA(size_t length, RNG& rng) {
	return std::string(length, DNA[rng()]);
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
	size_t lengthStart;
	size_t lengthEnd;
	size_t numTrials = 1;
	bool identical = false;
	bool different = false;
	bool repeats = false;
	boost::mt19937::result_type seed = 1;
	
	// Set up option parser
	util::options::Parser parser("",
								 "Calculates the newton polytope for the "
								 "polynomial representing all possible alignments "
								 "of the two sequences given as arguments.");
	parser.addStoreOpt('m', "match", "match variable", match, "SCORE");
	parser.addStoreOpt('x', "mismatch", "mismatch variable", mismatch, "SCORE");
	parser.addStoreOpt('s', "space", "space variable", space, "SCORE");
	parser.addStoreOpt('g', "gap", "gap variable", gap, "SCORE");
	parser.addStoreOpt(0, "matrix", "scoring matrix filename",
					   scoringMatrixFilename, "FILENAME");
	parser.addStoreOpt(0, "qhull-log",
					   "File for qhull warning/error messages",
					   qhullLogFile, "FILE");
	parser.addStoreOpt('t', "num-trials",
					   "Number of random trials for each sequence length",
					   numTrials, "NUM");
	parser.addStoreOpt('s', "seed",
					   "Seed for random number generator",
					   seed, "NUM");
	parser.addStoreConstOpt('i', "identical",
							"Make both sequences the same",
							identical, true);
	parser.addStoreConstOpt('r', "repeat",
							"Make sequences simple repeats",
							repeats, true);
	parser.addStoreConstOpt('d', "different",
							"Enforce that both sequences are different",
							different, true);
	parser.addStoreArg("lengthStart",
					   "Starting length of sequences to align",
					   lengthStart);
	parser.addStoreArg("lengthEnd",
					   "Ending length of sequences to align",
					   lengthEnd);
	parser.parse(argv, argv + argc);

	try {
		IntegerPolytope::setQhullLogFile(qhullLogFile);

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
				
		boost::mt19937 rng(seed);
		boost::uniform_int<> four(0,3);
		boost::variate_generator<boost::mt19937, boost::uniform_int<> >
			die(rng, four);

		boost::timer timer;
		
		for (size_t length = lengthStart; length <= lengthEnd; ++length) {
			size_t totalVertices = 0;
			double totalTime = 0;
			for (size_t trial = 0; trial < numTrials; ++trial) {
				std::string seq1, seq2;
				seq1 = (repeats ?
						simpleRepeatDNA(length, die) :
						randomDNA(length, die));
				if (identical) {
					seq2 = seq1;
				} else {
					do {
						seq2 = (repeats ?
								simpleRepeatDNA(length, die) :
								randomDNA(length, die));
					} while (different && length != 0 && seq1 == seq2);
				}
				timer.restart();
				IntegerPolytope p = calculatePolytope(seq1,
													  seq2,
													  matrix,
													  varMap[space],
													  varMap[gap],
													  numVars);

				totalVertices += p.getNumVertices();
				totalTime += timer.elapsed();
			}
			
			std::cout << length << '\t'
					  << static_cast<double>(totalVertices) / numTrials << '\t'
					  << totalTime / numTrials << std::endl;
		}

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
