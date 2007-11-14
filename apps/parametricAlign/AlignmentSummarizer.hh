#ifndef __ALIGNMENT_SUMMARIZER_HH__
#define __ALIGNMENT_SUMMARIZER_HH__

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "bio/alignment/ScoringMatrix.hh"
#include "bio/alignment/AlphabetScoringMatrix.hh"
#include "bio/alignment/AffineGapNWPairwiseSummarizer.hh"
using bio::alignment::ScoringMatrix;
using bio::alignment::AlphabetScoringMatrix;
using bio::alignment::AffineGapNWPairwiseSummarizer;
using bio::alignment::PairwiseAlignment;

template<typename Score>
class AlignmentSummarizer {
public:
	// Constructor from two sequencees to be aligned, a symbolic DNA
	// scoring matrix, and symbols for the space and gap scores
	AlignmentSummarizer(const std::string& seq1,
						const std::string& seq2,
						const ScoringMatrix<std::string>& varMatrix,
						const std::string spaceVar,
						const std::string gapVar,
						bool lexMin = false,
						bool linearMem = false,
						size_t maxMem =  (1 << 29));

	// Returns the number of parameters in the alignment model
	size_t getNumParams() const;

	// Returns the name of Ith parameter in the alignment model
	std::string getParam(size_t i) const;

	// Returns the lexicographically minimum optimal summary given
	// values for the alignment model parameters in VALUES
	std::vector<size_t> getSummary(const std::vector<Score>& values) const;

	std::vector<size_t> getSummary(const PairwiseAlignment& alignment) const;

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
	std::vector<std::string> symbolList;
	AlphabetScoringMatrix<size_t> indexMatrix;
	size_t spaceIndex;
	size_t gapIndex;
	bool isZeroParam;
	bool lexMin;
	bool linearMem;
	mutable AffineGapNWPairwiseSummarizer<Score> summarizer;
};

template<typename Score>
std::vector<size_t>
AlignmentSummarizer<Score>::
getSummary(const std::vector<Score>& values) const {
	summarizer.setScores(values);
	return summarizer.summarize(seq1, seq2, lexMin, linearMem);
}

template<typename Score>
std::vector<size_t>
AlignmentSummarizer<Score>::
getSummary(const PairwiseAlignment& alignment) const {
	std::vector<Score> values(getNumParams(), 0);
	summarizer.setScores(values);
	return summarizer.summarize(alignment);
}

template<typename Score>
const std::string AlignmentSummarizer<Score>::ZERO = "0";

template<typename Score>
AlignmentSummarizer<Score>::
AlignmentSummarizer(const std::string& seq1,
					const std::string& seq2,
					const ScoringMatrix<std::string>& varMatrix,
					const std::string spaceVar,
					const std::string gapVar,
					bool lexMin,
					bool linearMem,
					size_t maxMem)
	: seq1(seq1),
	  seq2(seq2),
	  indexMatrix(varMatrix.getAlphabet()),
	  lexMin(lexMin),
	  linearMem(linearMem),
	  summarizer(varMatrix.getAlphabet())
{
	makeSymbolList(varMatrix, spaceVar, gapVar);
	makeIndices(varMatrix, spaceVar, gapVar);
	summarizer.setIndices(indexMatrix, spaceIndex, gapIndex, isZeroParam);

	if (3 * (seq1.size() + 1) * (seq2.size() + 1) * sizeof(Score) > maxMem) {
		this->linearMem = true;
	}
}
	
template<typename Score>
size_t AlignmentSummarizer<Score>::getNumParams() const {
	return symbolList.size();
}

template<typename Score>
std::string AlignmentSummarizer<Score>::getParam(size_t i) const {
	return symbolList[i];
}

template<typename Score>
void
AlignmentSummarizer<Score>::
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
AlignmentSummarizer<Score>::
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

#endif // __ALIGNMENT_SUMMARIZER_HH__
