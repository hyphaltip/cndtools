#include "bio/alignment/AlphabetScoringMatrix.hh"

#include "LikelihoodCalculator.hh"
#include "RecombinationScorer.hh"
#include "LogLikelihoodSemiRing.hh"

const bio::alphabet::Nucleotide
LikelihoodCalculator::GAPPED_DNA("ACGTNMRWSYKVHDB-",
								 "TGCANKYWSRMBDHV-");

LogDouble LikelihoodCalculator::likelihood(const MultipleAlignment& ma,
										   size_t recombinant,
										   double matchProb,
										   double norecombinationProb) {
	using bio::alignment::AlphabetScoringMatrix;
	
	typedef LogLikelihoodSemiRing SemiRing;
	typedef SemiRing::Element Element;
	SemiRing semiRing;

	const int N = ma.getNumSeqs() - 1;

	Element match(matchProb);
	Element mismatch((1.0 - matchProb) / 4);
	Element norecombination(norecombinationProb);
	Element recombination((1.0 - norecombinationProb) / (N - 1));
	Element initial(1.0 / N);

	std::vector<Element> initialScores(N + 1, initial);
	AlphabetScoringMatrix<Element> emissionScores(GAPPED_DNA, match, mismatch);
	util::Matrix<Element> transitionScores(N + 1, N + 1);
	for (int i = 0; i <= N; ++i) {
		for (int j = 0; j <= N; ++j) {
			transitionScores(i, j) = i == j ? norecombination : recombination;
		}
	}

	RecombinationScorer<SemiRing>
		scorer(semiRing, initialScores, emissionScores, transitionScores);

	return scorer.score(ma, recombinant);
}
