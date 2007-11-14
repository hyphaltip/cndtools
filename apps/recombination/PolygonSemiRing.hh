#include "Polygon.hh"

typedef Polygon<int> IntegerPolygon;
typedef Vector<int> IntegerVector;

class PolygonSemiRing {
public:
	typedef IntegerPolygon Element;

	PolygonSemiRing() {}

	Element getZero() const { return IntegerPolygon(); }
	Element getMultiplicativeIdentity() const { return IntegerVector(2); }
};
