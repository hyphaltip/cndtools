#include "math/LogDouble.hh"

class LogLikelihoodSemiRing {
public:
	typedef math::LogDouble Element;

	LogLikelihoodSemiRing() {}

	Element getZero() const { return 0.0; }
	Element getMultiplicativeIdentity() const { return 1.0; }
};
