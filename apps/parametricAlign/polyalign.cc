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

#include <stdexcept>

#include "polyalign.hh"

#include "bio/alignment/NeedlemanWunschPairwiseScorer.hh"
#include "bio/alignment/AffineGapNWPairwiseScorer.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "bio/alignment/DNAScoringMatrix.hh"
#include "util/string.hh"
using namespace bio::alignment;
using namespace polytope;
using util::string::toString;
using bio::alphabet::DNA;

class PolytopeSemiRing {
public:
	typedef Polytope<int> Element;

	PolytopeSemiRing(size_t dims) : dims(dims) {}

	Element getZero() const { return Polytope<int>(); }
	Element getMultiplicativeIdentity() const { return Vector<int>(dims); }

private:
	size_t dims;
};

IntegerPolytope calculatePolytope(const std::string& s1,
								  const std::string& s2,
								  PolytopeMatrix& matrix,
								  IntegerPolytope& space,
								  IntegerPolytope& gap,
								  size_t numVars) {
	typedef PolytopeSemiRing SemiRing;
	SemiRing semiRing(numVars);
	AffineGapNWPairwiseScorer<SemiRing>
		scorer(semiRing, matrix, space, gap);
	return scorer.score(s1, s2);
}	

static void
calcMiddleRowPolytopesGlobal4(const std::string& seq1,
							  const std::string& seq2,
							  PolytopeList& forwardMatch,
							  PolytopeList& forwardGap1,
							  PolytopeList& forwardGap2,
							  PolytopeList& backwardMatch,
							  PolytopeList& backwardGap1,
							  PolytopeList& backwardGap2) {
	enum PARAMS {MISMATCH, SPACE, GAP, NUM_PARAMS};
	Polytope<int>
		match = Vector<int>(NUM_PARAMS),
		mismatch = Vector<int>::unit(NUM_PARAMS, MISMATCH),
		space = Vector<int>::unit(NUM_PARAMS, SPACE),
		gap = Vector<int>::unit(NUM_PARAMS, GAP);

	typedef PolytopeSemiRing SemiRing;
	SemiRing semiRing(NUM_PARAMS);
	DNAScoringMatrix< Polytope<int> > matrix(match, mismatch);
	AffineGapNWPairwiseScorer<SemiRing>
		scorer(semiRing, matrix, space, gap);
	size_t seq1MiddlePos = seq1.size() / 2;
	std::string seq1Prefix = seq1.substr(0, seq1MiddlePos);
	std::string seq1Suffix = seq1.substr(seq1MiddlePos);
	scorer.scoreLastRow(seq1Prefix,
						seq2,
						forwardMatch,
						forwardGap1,
						forwardGap2);
	scorer.scoreFirstRow(seq1Suffix,
						 seq2,
						 backwardMatch,
						 backwardGap1,
						 backwardGap2);
}

static void
calcMiddleRowPolytopesGlobal5(const std::string& seq1,
							  const std::string& seq2,
							  PolytopeList& forwardMatch,
							  PolytopeList& forwardGap1,
							  PolytopeList& forwardGap2,
							  PolytopeList& backwardMatch,
							  PolytopeList& backwardGap1,
							  PolytopeList& backwardGap2) {
	enum PARAMS {TRANSITION, TRANSVERSION, SPACE, GAP, NUM_PARAMS};
	Polytope<int>
		match = Vector<int>(NUM_PARAMS),
		transition = Vector<int>::unit(NUM_PARAMS, TRANSITION),
		transversion = Vector<int>::unit(NUM_PARAMS, TRANSVERSION),
		space = Vector<int>::unit(NUM_PARAMS, SPACE),
		gap = Vector<int>::unit(NUM_PARAMS, GAP);

	typedef PolytopeSemiRing SemiRing;
	SemiRing semiRing(NUM_PARAMS);
	DNAScoringMatrix< Polytope<int> > matrix(match, transition, transversion);
	AffineGapNWPairwiseScorer<SemiRing>
		scorer(semiRing, matrix, space, gap);
	size_t seq1MiddlePos = seq1.size() / 2;
	std::string seq1Prefix = seq1.substr(0, seq1MiddlePos);
	std::string seq1Suffix = seq1.substr(seq1MiddlePos);
	scorer.scoreLastRow(seq1Prefix,
						seq2,
						forwardMatch,
						forwardGap1,
						forwardGap2);
	scorer.scoreFirstRow(seq1Suffix,
						 seq2,
						 backwardMatch,
						 backwardGap1,
						 backwardGap2);
}

void
calcMiddleRowPolytopes(const std::string& seq1,
					   const std::string& seq2,
					   size_t num_params,
					   PolytopeList& forwardMatch,
					   PolytopeList& forwardGap1,
					   PolytopeList& forwardGap2,
					   PolytopeList& backwardMatch,
					   PolytopeList& backwardGap1,
					   PolytopeList& backwardGap2) {
	if (num_params == 4) {
		return calcMiddleRowPolytopesGlobal4(seq1,
											 seq2,
											 forwardMatch,
											 forwardGap1,
											 forwardGap2,
											 backwardMatch,
											 backwardGap1,
											 backwardGap2);
	} else if (num_params == 5) {
		return calcMiddleRowPolytopesGlobal5(seq1,
											 seq2,
											 forwardMatch,
											 forwardGap1,
											 forwardGap2,
											 backwardMatch,
											 backwardGap1,
											 backwardGap2);
	} else {
		throw std::runtime_error("Bad number of parameters: " +
								 toString(num_params));
	}
}

const std::string ZERO = "0";

VariableSet
makeVariableSet(VariableMatrix& matrix,
				Variable& space,
				Variable& gap) {
	// Collect all of the variable names
	VariableSet varSet;
	varSet.insert(space);
	varSet.insert(gap);
	for (size_t i = 0; i < DNA.getSize(); ++i) {
		for (size_t j = 0; j < DNA.getSize(); ++j) {
			std::string entry = matrix.getScore(i, j);
			varSet.insert(entry);
		}
	}
	return varSet;
}

size_t
getNumVars(const VariableSet& varSet) {
	bool zeroIncluded = (varSet.find(ZERO) != varSet.end());
	return (zeroIncluded ? varSet.size() - 1 : varSet.size());
}

VariablePolytopeMap
makeVariablePolytopeMap(VariableSet& varSet) {
	VariablePolytopeMap varMap;
	// Determine the number of variables and if some parameters are zero
	bool zeroIncluded = (varSet.find(ZERO) != varSet.end());
	size_t numVars = (zeroIncluded ? varSet.size() - 1 : varSet.size());

	// Create variable polytopes
	size_t i = 0;
	typedef VariableSet::const_iterator VarIt;
	for (VarIt v = varSet.begin(); v != varSet.end(); ++v) {
		if (*v == ZERO) { continue; }
		varMap[*v] = Vector<int>::unit(numVars, i);
		++i;
	}

	// Add zero if necessary
	if (zeroIncluded) {
		varMap[ZERO] = Vector<int>(numVars);
	}

	return varMap;
}

PolytopeMatrix
makePolytopeMatrix(VariableMatrix& varMatrix,
				   VariablePolytopeMap& varMap) {
	PolytopeMatrix matrix;
	for (size_t i = 0; i < DNA.getSize(); ++i) {
		for (size_t j = 0; j < DNA.getSize(); ++j) {
			matrix.setScore(i, j, varMap[varMatrix.getScore(i, j)]);
		}
	}
	return matrix;
}

VariableDimensionMap
makeVariableDimensionMap(VariableSet& varSet) {
	VariableDimensionMap varDimMap;
	size_t i = 1;
	typedef VariableSet::const_iterator VarIt;
	for (VarIt v = varSet.begin(); v != varSet.end(); ++v) {
		if (*v == ZERO) {
			varDimMap[*v] = 0;
		} else {
			varDimMap[*v] = i;
			++i;
		}
	}
	return varDimMap;
}

void
setParameters(IntegerMatrix& scoreMatrix, int& spaceScore, int& gapScore,
			  VariableMatrix& varMatrix, Variable& space, Variable& gap,
			  VariableSet& varSet,
			  IntegerPolytope::Ray& paramRay) {
	VariableDimensionMap varDimMap = makeVariableDimensionMap(varSet);

	// Set matrix scores
	for (size_t i = 0; i < DNA.getSize(); ++i) {
		for (size_t j = 0; j < DNA.getSize(); ++j) {
			size_t dim = varDimMap[varMatrix.getScore(i, j)];
			scoreMatrix.setScore(i, j, paramRay[dim]);
		}
	}

	// Set gap and space scores
	spaceScore = paramRay[varDimMap[space]];
	gapScore = paramRay[varDimMap[gap]];
}
