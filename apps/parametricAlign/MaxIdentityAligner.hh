#ifndef __MAXIDENTITYALIGNER_HH__
#define __MAXIDENTITYALIGNER_HH__

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "util/stl.hh"
#include "bio/alignment/ScoringMatrix.hh"
#include "bio/alignment/AlphabetScoringMatrix.hh"
#include "bio/alignment/AffineGapNWPairwiseSummarizer.hh"
#include "bio/alignment/AffineGapNWPairwiseAligner.hh"
#include "bio/alignment/AffineGapNWFastPairwiseAligner.hh"
#include "MaxIdentityPairwiseAligner.hh"
using bio::alignment::ScoringMatrix;
using bio::alignment::AlphabetScoringMatrix;
using bio::alignment::VarMappingScoringMatrix;
using bio::alignment::AffineGapNWPairwiseAligner;
using bio::alignment::AffineGapNWFastPairwiseAligner;
using bio::alignment::AffineGapNWPairwiseSummarizer;
using bio::alignment::PairwiseAlignment;

template<typename Score>
class MaxIdentityAligner {
public:
	// Constructor from two sequencees to be aligned, a symbolic DNA
	// scoring matrix, and symbols for the space and gap scores
	MaxIdentityAligner(const std::string& seq1,
					   const std::string& seq2,
					   size_t start,
					   size_t end,
					   const ScoringMatrix<std::string>& varMatrix,
					   const std::string spaceVar,
					   const std::string gapVar);
	
	// Returns the number of parameters in the alignment model
	size_t getNumParams() const;

	// Returns the name of Ith parameter in the aligment model
	std::string getParam(size_t i) const;

	// Returns the summary of the subalignment with the maximum number
	// of matches
	PairwiseAlignment getAlignment(const std::vector<Score>& values) const;

private:
	void makeSymbolList(const ScoringMatrix<std::string>& varMatrix,
						const std::string spaceVar,
						const std::string gapVar);
	void makeIndices(const ScoringMatrix<std::string>& varMatrix,
					 const std::string spaceVar,
					 const std::string gapVar);

	static const std::string ZERO;

	std::string seq1;
	std::string seq2;
	size_t start;
	size_t end;
	std::vector<std::string> symbolList;
	AlphabetScoringMatrix<size_t> indexMatrix;
	size_t spaceIndex;
	size_t gapIndex;
	mutable MaxIdentityPairwiseAligner<Score> aligner;
	mutable AffineGapNWPairwiseSummarizer<Score> summarizer;
	bool isZeroParam;
};

template<typename Score>
PairwiseAlignment
MaxIdentityAligner<Score>::
getAlignment(const std::vector<Score>& values) const {
	std::vector<Score> eltScores(values);
	if (isZeroParam) {
		eltScores.push_back(0);
	}
// 	summarizer.setScores(values);
	
	VarMappingScoringMatrix<Score> matrix(indexMatrix, eltScores);
	aligner.setScores(matrix, eltScores[spaceIndex], eltScores[gapIndex]);
	PairwiseAlignment alignment = aligner.align(seq1, seq2, start, end);

	PairwiseAlignment slice = alignment.slice(0, start, end);
	return slice;
}

template<typename Score>
const std::string MaxIdentityAligner<Score>::ZERO = "0";

template<typename Score>
MaxIdentityAligner<Score>::
MaxIdentityAligner(const std::string& seq1,
				   const std::string& seq2,
				   size_t start,
				   size_t end,
				   const ScoringMatrix<std::string>& varMatrix,
				   const std::string spaceVar,
				   const std::string gapVar)
	: seq1(seq1),
	  seq2(seq2),
	  start(start),
	  end(end),
	  indexMatrix(varMatrix.getAlphabet()),
	  aligner(varMatrix.getAlphabet()),
	  summarizer(varMatrix.getAlphabet())
{
	makeSymbolList(varMatrix, spaceVar, gapVar);
	makeIndices(varMatrix, spaceVar, gapVar);
	summarizer.setIndices(indexMatrix, spaceIndex, gapIndex, isZeroParam);
}
	
template<typename Score>
size_t MaxIdentityAligner<Score>::getNumParams() const {
	return symbolList.size();
}

template<typename Score>
std::string MaxIdentityAligner<Score>::getParam(size_t i) const {
	return symbolList[i];
}

template<typename Score>
void
MaxIdentityAligner<Score>::
makeSymbolList(const ScoringMatrix<std::string>& varMatrix,
			   const std::string spaceVar,
			   const std::string gapVar) {
	size_t alphabetSize = varMatrix.getAlphabet().getSize();
	// Add matrix symbols to list
	for (size_t i = 0; i < alphabetSize; ++i) {
		for (size_t j = 0; j < alphabetSize; ++j) {
			std::string entry = varMatrix.getScore(i, j);
			symbolList.push_back(entry);
		}
	}
	// Add space and gap symbols to list
	symbolList.push_back(spaceVar);
	symbolList.push_back(gapVar);

	// Make ordered unique list of symbols
	std::sort(symbolList.begin(), symbolList.end());
	symbolList.erase(std::unique(symbolList.begin(), symbolList.end()),
					 symbolList.end());

	// Check if ZERO is a symbol, if so remove it and note its presence
	std::vector<std::string>::iterator zeroPos =
		std::find(symbolList.begin(), symbolList.end(), "0");
	if (zeroPos != symbolList.end()) {
		isZeroParam = true;
		symbolList.erase(zeroPos);
	} else {
		isZeroParam = false;
	}
}

template<typename Score>
void
MaxIdentityAligner<Score>::
makeIndices(const ScoringMatrix<std::string>& varMatrix,
			const std::string spaceVar,
			const std::string gapVar) {
	// Make a map between symbols and indices
	std::map<std::string, size_t> symIndexMap;
	for (size_t i = 0; i < symbolList.size(); ++i) {
		symIndexMap[symbolList[i]] = i;
	}

	// Add an index for ZERO, if necessary
	if (isZeroParam) {
		symIndexMap[ZERO] = symbolList.size();
	}

	// Set entries in index matrix
	size_t alphabetSize = varMatrix.getAlphabet().getSize();
	for (size_t i = 0; i < alphabetSize; ++i) {
		for (size_t j = 0; j < alphabetSize; ++j) {
			std::string s = varMatrix.getScore(i, j);
			indexMatrix.setScore(i, j, symIndexMap[s]);
		}
	}

	// Set space and gap indices
	spaceIndex = symIndexMap[spaceVar];
	gapIndex = symIndexMap[gapVar];
}

#endif // __MAXIDENTITYALIGNER_HH__
