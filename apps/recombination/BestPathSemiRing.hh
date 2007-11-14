#include <vector>
#include <algorithm>

#include "boost/shared_ptr.hpp"
using boost::shared_ptr;

struct Path;
typedef shared_ptr<Path> PathPtr;

struct Path {
	size_t index;
	PathPtr last;

	Path(size_t index = 0, PathPtr last = PathPtr())
		: index(index), last(last) {}

	Path(PathPtr head, PathPtr tail) {
		index = head->index;
		if (head->last.get() == NULL) {
			last = tail;
		} else {
			last = PathPtr(new Path(head->last, tail));
		}
	}
		
};

void get_path_indices(PathPtr path, std::vector<size_t>& indices) {
	while (path != NULL) {
		indices.push_back(path->index);
		path = path->last;
	}
	std::reverse(indices.begin(), indices.end());
}

template<typename Score>
struct ScoredPath {
	Score score;
	PathPtr path;

	ScoredPath(Score score = Score(), PathPtr path = PathPtr())
		: score(score), path(path) {}

	ScoredPath& operator+=(const ScoredPath& x) {
		if (score < x.score) { *this = x; }
		return *this;
	}

	ScoredPath& operator*=(const ScoredPath& x) {
		score *= x.score;
		if (x.path.get() != NULL) {
			if (path.get() != NULL) {
				path = PathPtr(new Path(x.path, path));
			} else {
				path = x.path;
			}
		}
		return *this;
	}

	ScoredPath operator*(const ScoredPath& x) {
		return ScoredPath(*this) *= x;
	}

	ScoredPath operator+(const ScoredPath& x) {
		return ScoredPath(*this) += x;
	}
};

template<typename ScoringSemiRing>
class BestPathSemiRing {
public:
	typedef ScoredPath<typename ScoringSemiRing::Element> Element;

	BestPathSemiRing(const ScoringSemiRing& scoringSemiRing)
		: scoringSemiRing(scoringSemiRing) {
	}
	
	Element getZero() const {
		return Element(scoringSemiRing.getZero());
	}
	
	Element getMultiplicativeIdentity() const {
		return Element(scoringSemiRing.getMultiplicativeIdentity());
	}

private:
	ScoringSemiRing scoringSemiRing;
};
