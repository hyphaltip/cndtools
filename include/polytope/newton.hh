/* Copyright (c) 2006
   Colin Dewey (University of Wisconsin-Madison)
   cdewey@biostat.wisc.edu
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __NEWTON_HH__
#define __NEWTON_HH__

#include <ostream>

#include <Poly.h>
#include <Rational.h>
#include <ListMatrix.h>
#include <Matrix.h>
#include <Vector.h>
#include <Set.h>
#include <Array.h>

namespace polytope {

	using namespace polymake;

	class Newton {		
	protected:
		
		typedef ListMatrix< Vector<Rational> > VerticesType;
		mutable VerticesType vertices;
		mutable bool reduced;
		typedef Matrix<Rational> FacetsType;
		mutable FacetsType facets;
		mutable bool facetsComputed;
		mutable int lastReducedSize;

		static const int POINT_FACTOR;
		
		void reduce() const;
		void reduceIfBig() const;
		void computeFacets() const;
	
	public:
		Newton(int dims = 0);
		Newton(const Newton& other);

		const VerticesType& getVertices() const;
		const size_t getNumVertices() const;

		const FacetsType& getFacets() const;
		const size_t getNumFacets() const;


		Array< Set<int> > getFacetsThruVertices();
		Array<Poly> getNormalFan();
		Array<Vector<Integer> > getFacetNormals();
		
		Newton& operator+=(const Newton& other);
		Newton& operator+=(Rational r);
		Newton& operator*=(const Newton& other);
		Newton& operator*=(Rational r);
		Newton& operator^=(Rational r);

		Newton operator+(const Newton& other) const;
		Newton operator+(Rational r) const;		
		Newton operator*(const Newton& other) const;
		Newton operator*(Rational r) const;		
		Newton operator^(Rational r) const;

		friend std::ostream& operator<<(std::ostream& strm, const Newton& p);
										
	};

	class Variable : public Newton {
	public:
		Variable(int varIndex, int numVars);
	};

	class Constant : public Newton {
	public:
		Constant(Rational r, int numVars);
	};

	inline Newton::Newton(int dims)
		: vertices(0, dims + 1),
		  reduced(true),
		  facetsComputed(false),
		  lastReducedSize(0)
	{}
	
	inline Newton::Newton(const Newton& other)
		: vertices(other.vertices),
		  reduced(other.reduced),
		  facetsComputed(other.facetsComputed),
		  lastReducedSize(other.lastReducedSize)
	{}
	
	inline Newton Newton::operator+(const Newton& other) const {
		Newton result(*this);
		result += other;
		return result;
	}

	inline Newton Newton::operator+(Rational r) const {
		Newton result(*this);
		result += r;
		return result;
	}
	
	inline Newton Newton::operator*(const Newton& other) const {
		Newton result(*this);
		result *= other;
		return result;
	}

	inline Newton Newton::operator*(Rational r) const {
		Newton result(*this);
		result *= r;
		return result;
	}
	
 	inline Newton Newton::operator^(Rational r) const {
		Newton result(*this);
		result ^= r;
		return result;
	}
	
	inline std::ostream& operator<<(std::ostream& strm, const Newton& p) {
		wrap(strm) << p.getVertices();
		return strm;
	}
};

#endif // __NEWTON_HH__
