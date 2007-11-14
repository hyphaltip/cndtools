#ifndef __LIKELIHOODCALCULATOR_HH__
#define __LIKELIHOODCALCULATOR_HH__

#include "bio/alphabet/Nucleotide.hh"
#include "bio/alignment/MultipleAlignment.hh"
#include "math/LogDouble.hh"
using math::LogDouble;
using bio::alignment::MultipleAlignment;

class LikelihoodCalculator {
public:
	
	LikelihoodCalculator() {}
	
	static LogDouble likelihood(const MultipleAlignment& ma,
								size_t recombinant,
								double matchProb,
								double norecombinationProb);

private:
	static const bio::alphabet::Nucleotide GAPPED_DNA;
};

#endif // __LIKELIHOODCALCULATOR_HH__
