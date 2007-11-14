#ifndef __BESTSTATESETLISTER_HH__
#define __BESTSTATESETLISTER_HH__

#include "bio/alignment/AlphabetScoringMatrix.hh"
#include "bio/alphabet/Nucleotide.hh"
#include "math/MaxPlus.hh"
using bio::alignment::AlphabetScoringMatrix;

#include "RecombinationScorer.hh"
#include "BestStateSetListSemiRing.hh"

class BestStateSetLister {
public:
	typedef int Score;

	BestStateSetLister(Score matchScore = Score(),
					   Score mismatchScore = Score(),
					   Score recombinationScore = Score(),
					   Score norecombinationScore = Score());

	void getBestStateSetList(const MultipleAlignment& ma,
							 size_t recombinant,
							 StateSetList& stateSetList);

	void setScores(Score matchScore,
				   Score mismatchScore,
				   Score recombinationScore,
				   Score norecombinationScore);

private:
	static const bio::alphabet::Nucleotide GAPPED_DNA;
	
	Score matchScore;
	Score mismatchScore;
	Score recombinationScore;
	Score norecombinationScore;
};

const bio::alphabet::Nucleotide
BestStateSetLister::GAPPED_DNA("ACGTNMRWSYKVHDB-",
							   "TGCANKYWSRMBDHV-");

BestStateSetLister::
BestStateSetLister(Score matchScore, Score mismatchScore,
				   Score recombinationScore, Score norecombinationScore)
	: matchScore(matchScore),
	  mismatchScore(mismatchScore),
	  recombinationScore(recombinationScore),
	  norecombinationScore(norecombinationScore) {
}

void
BestStateSetLister::
setScores(Score matchScore,
		  Score mismatchScore,
		  Score recombinationScore,
		  Score norecombinationScore) {
	this->matchScore = matchScore;
	this->mismatchScore = mismatchScore;
	this->recombinationScore = recombinationScore;
	this->norecombinationScore = norecombinationScore;
}

void
BestStateSetLister::
getBestStateSetList(const MultipleAlignment& ma,
					size_t recombinant,
					StateSetList& stateSetList) {
	typedef BestStateSetListSemiRing< math::MaxPlus<Score> > SemiRing;
	typedef SemiRing::Element Element;
	math::MaxPlus<Score> scoringSemiRing;
	SemiRing semiRing(scoringSemiRing, ma.getNumSeqs());
	StateSet emptySet(ma.getNumSeqs());
	StateSetList emptySetList(1, emptySet);
	Element match(matchScore, emptySetList);
	Element mismatch(mismatchScore, emptySetList);
	std::vector<Element> initialScores(ma.getNumSeqs(),
									   semiRing.getMultiplicativeIdentity());
	AlphabetScoringMatrix<Element> emissionScores(GAPPED_DNA, match, mismatch);
	util::Matrix<Element> transitionScores(ma.getNumSeqs(),
										   ma.getNumSeqs());
	for (size_t i = 0; i < ma.getNumSeqs(); ++i) {
		for (size_t j = 0; j < ma.getNumSeqs(); ++j) {
			StateSet s(ma.getNumSeqs());
			s.set(j);
			StateSetList stateSetList(1, s);
			Score score = (i == j ? norecombinationScore : recombinationScore);
			transitionScores(i, j) = Element(score, stateSetList);
		}
	}

	RecombinationScorer<SemiRing>
		scorer(semiRing, initialScores, emissionScores, transitionScores);

	Element best = scorer.score(ma, recombinant);
	stateSetList = best.stateSetList;
}

#endif // __BESTSTATESETLISTER_HH__
