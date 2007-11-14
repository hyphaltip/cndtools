#ifndef __MAXIDENTITYPAIRWISEALIGNER_HH__
#define __MAXIDENTITYPAIRWISEALIGNER_HH__

#include "bio/alignment/AlphabetScoringMatrix.hh"
#include "bio/alphabet/Alphabet.hh"
#include "bio/alignment/PairwiseAlignment.hh"
using namespace bio;
using namespace bio::alignment;

template<typename NumType>
struct Score {
	NumType val;
	size_t identities;

	Score(NumType val = 0, size_t identities = 0) :
		val(val), identities(identities) {
	}

	Score& operator=(NumType val) {
		this->val = val;
		this->identities = 0;
		return *this;
	}

	Score& operator+=(const Score& s) {
		val += s.val;
		identities += s.identities;
		return *this;
	}

	Score operator+(const Score& s) const {
		Score t(*this);
		return t += s;
	}

	bool operator<(const Score& s) const {
		return val < s.val or val == s.val and identities < s.identities;
	}

};

template<typename NumType>
std::ostream& operator<<(std::ostream& stream, const Score<NumType>& s) {
	return stream << '(' << s.val << ',' << s.identities << ')';
}

template<typename NumType>	
class MaxIdentityPairwiseAligner {
public:
	MaxIdentityPairwiseAligner(const alphabet::Alphabet& alphabet);
		
	MaxIdentityPairwiseAligner(const ScoringMatrix<NumType>& matrix,
							   const NumType& space,
							   const NumType& gap);
		
	PairwiseAlignment align(const std::string& seq1,
							const std::string& seq2,
							size_t start,
							size_t end) const;
	
	void setScores(const ScoringMatrix<NumType>& matrix,
				   const NumType& space,
				   const NumType& gap);
		
private:
	static PairwiseAlignment
	statesToAlignment(const std::vector<size_t>& indices,
					  const std::string& seq1,
					  const std::string& seq2);

	static size_t argmax(Score<NumType> h,
						 Score<NumType> d,
						 Score<NumType> i);

	static const size_t H;
	static const size_t D;
	static const size_t I;

	AlphabetScoringMatrix<NumType> matrix;
	Score<NumType> space;
	Score<NumType> gap;
	mutable util::Matrix< Score<NumType> > hState;
	mutable util::Matrix< Score<NumType> > dState;
	mutable util::Matrix< Score<NumType> > iState;
};

template<typename NumType>	
const size_t MaxIdentityPairwiseAligner<NumType>::H = 0;
template<typename NumType>	
const size_t MaxIdentityPairwiseAligner<NumType>::D = 1;
template<typename NumType>	
const size_t MaxIdentityPairwiseAligner<NumType>::I = 2;

template<typename NumType>
MaxIdentityPairwiseAligner<NumType>::
MaxIdentityPairwiseAligner(const alphabet::Alphabet& alphabet) 
	: matrix(alphabet),
	  space(),
	  gap() {
}
	
template<typename NumType>
MaxIdentityPairwiseAligner<NumType>::
MaxIdentityPairwiseAligner(const ScoringMatrix<NumType>& matrix,
							   const NumType& space,
							   const NumType& gap)
	: matrix(matrix),
	  space(space),
	  gap(gap) {
}

template<typename NumType>
void
MaxIdentityPairwiseAligner<NumType>::
setScores(const ScoringMatrix<NumType>& matrix,
		  const NumType& space,
		  const NumType& gap) {
	this->matrix = matrix;
	this->space = Score<NumType>(space);
	this->gap = Score<NumType>(gap);
}
	
template<typename NumType>
size_t
MaxIdentityPairwiseAligner<NumType>::
argmax(Score<NumType> h, Score<NumType> d, Score<NumType> i) {
	if (h < d and i < d) {
		return D;
	} else if (h < i) {
		return I;
	} else {
		return H;
	}
}

template<typename NumType>
PairwiseAlignment
MaxIdentityPairwiseAligner<NumType>::
statesToAlignment(const std::vector<size_t>& indices,
				  const std::string& seq1,
				  const std::string& seq2) {
	PairwiseAlignment alignment;
	size_t pos1 = 0, pos2 = 0;
	for (size_t i = 0; i < indices.size(); ++i) {
		if (indices[i] == H) {
			alignment.seq1 += seq1[pos1];
			alignment.seq2 += seq2[pos2];
			++pos1;
			++pos2;
		} else if (indices[i] == I) {
			alignment.seq1 += '-';
			alignment.seq2 += seq2[pos2];
			++pos2;
		} else {
			alignment.seq1 += seq1[pos1];
			alignment.seq2 += '-';
			++pos1;
		}
	}

	return alignment;
}
	
template<typename NumType>
PairwiseAlignment
MaxIdentityPairwiseAligner<NumType>::
align(const std::string& seq1,
	  const std::string& seq2,
	  size_t start,
	  size_t end) const {
	using std::max;

	std::string eSeq1(matrix.getAlphabet().encode(seq1));
	std::string eSeq2(matrix.getAlphabet().encode(seq2));
		
	size_t n = seq1.size();
	size_t m = seq2.size();

	hState.resize(n + 1, m + 1);
	dState.resize(n + 1, m + 1);
	iState.resize(n + 1, m + 1);
		
	hState(n, m) = dState(n, m) = iState(n, m) = 0;

	// Calculate last row
	for (size_t j = 0; j < m; ++j) {
		size_t s = m - j; // source column
		size_t t = s - 1; // target column
		iState(n, t) = iState(n, s) + space;
		hState(n, t) = dState(n, t) = iState(n, t) + gap;
	}

	// Calculate all rows
	for (size_t i = 1; i <= n; ++i) {
		size_t tRow = n - i; // target row
		size_t sRow = tRow + 1; // source row
		for (size_t j = 0; j <= m; ++j) {
			size_t tCol = m - j; // target column
			size_t sCol = tCol + 1; // source column

			Score<NumType> dVal = dState(sRow, tCol) + space;
			Score<NumType> hVal = dVal + gap;
			Score<NumType> iVal = hVal;

			if (tCol < m) {
				Score<NumType> diag = hState(sRow, sCol)
					+ Score<NumType>(matrix.getScore(eSeq1[tRow], eSeq2[tCol]),
									 (eSeq1[tRow] == eSeq2[tCol] and
									  tRow >= start and
									  tRow < end));
				hVal = iVal = max(diag, hVal);
				dVal = max(diag, dVal);
					
				Score<NumType> horiz = iState(tRow, sCol) + space;
				hVal = max(horiz + gap, hVal);
				dVal = max(horiz + gap, dVal);
				iVal = max(horiz, iVal);
			}

			hState(tRow, tCol) = hVal;
			dState(tRow, tCol) = dVal;
			iState(tRow, tCol) = iVal;				
		}
	}

	//std::cerr << hState(0, 0) << '\n';
	
	std::vector<size_t> traceback;

	size_t i = 0, j = 0;
		
	size_t currState = H;

	while (i < n or j < m) {
		if (i == n) {
			currState = I;
		} else if (j == m) {
			currState = D;
		} else {
			Score<NumType> hScore = hState(i + 1, j + 1)
				+ Score<NumType>(matrix.getScore(eSeq1[i], eSeq2[j]),
								 (eSeq1[i] == eSeq2[j] and 
								  i >= start and
								  i < end));
			Score<NumType> dScore = dState(i + 1, j) + space;
			Score<NumType> iScore = iState(i, j + 1) + space;

			if (currState != D) { dScore += gap; }
			if (currState != I) { iScore += gap; }

			currState = argmax(hScore, dScore, iScore);
		}

		switch (currState) {
		case H:
			++i;
			++j;
			break;
		case D:
			++i;
			break;
		case I:
			++j;
			break;
		}

		traceback.push_back(currState);
	}
		
	return statesToAlignment(traceback, seq1, seq2);
}


#endif // __MAXIDENTITYPAIRWISEALIGNER_HH__
