#include <vector>
#include <algorithm>

#include "boost/dynamic_bitset.hpp"
using boost::dynamic_bitset;

typedef dynamic_bitset<> StateSet;
typedef std::vector<StateSet> StateSetList;

template<typename Score>
struct ScoredStateSetList {
	Score score;
	StateSetList stateSetList;

	ScoredStateSetList(Score score = Score(),
					   const StateSetList& stateSetList = StateSetList())
		: score(score), stateSetList(stateSetList) {}
	
	ScoredStateSetList& operator+=(const ScoredStateSetList& x) {
		if (score < x.score) {
			*this = x;
		} else if (score == x.score) {
			StateSetList unionList;
			std::set_union(stateSetList.begin(), stateSetList.end(),
						   x.stateSetList.begin(), x.stateSetList.end(),
						   std::back_inserter(unionList));
			stateSetList.swap(unionList);
		}
		return *this;
	}

	ScoredStateSetList& operator*=(const ScoredStateSetList& x) {
		score *= x.score;
		if (x.stateSetList.size() == 1) {
			StateSet multiplicand(x.stateSetList.front());
			if (multiplicand.none()) { return *this; }
			for (size_t i = 0; i < stateSetList.size(); ++i) {
				stateSetList[i] |= multiplicand;
			}
		} else {
			StateSetList product;
			for (size_t i = 0; i < stateSetList.size(); ++i) {
				for (size_t j = 0; j < x.stateSetList.size(); ++j) {
					product.push_back(stateSetList[i] |
									  x.stateSetList[j]);
				}
			}
			stateSetList.swap(product);
		}
		std::sort(stateSetList.begin(), stateSetList.end());
		stateSetList.erase(std::unique(stateSetList.begin(), stateSetList.end()),
						   stateSetList.end());
		return *this;
	}

	ScoredStateSetList operator*(const ScoredStateSetList& x) {
		return ScoredStateSetList(*this) *= x;
	}

	ScoredStateSetList operator+(const ScoredStateSetList& x) {
		return ScoredStateSetList(*this) += x;
	}
};

template<typename ScoringSemiRing>
class BestStateSetListSemiRing {
public:
	typedef ScoredStateSetList<typename ScoringSemiRing::Element> Element;

	BestStateSetListSemiRing(const ScoringSemiRing& scoringSemiRing,
							 size_t numStates)
		: scoringSemiRing(scoringSemiRing), numStates(numStates) {
	}
	
	Element getZero() const {
		return Element(scoringSemiRing.getZero());
	}
	
	Element getMultiplicativeIdentity() const {
		StateSet noStates(numStates);
		StateSetList stateSetList(1, noStates);
		return Element(scoringSemiRing.getMultiplicativeIdentity(),
								  stateSetList);
	}

private:
	ScoringSemiRing scoringSemiRing;
	size_t numStates;
};
