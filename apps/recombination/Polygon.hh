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

#ifndef __POLYGON_HH__
#define __POLYGON_HH__

#include <vector>
#include <algorithm>
#include <iosfwd>

#include "Vector.hh"

template<typename T>
class Polygon {
public:
	typedef typename std::vector< Vector<T> > VectorList;

	Polygon();
	Polygon(const Vector<T>& v);
		
	Polygon& operator+=(const Polygon& other);
	Polygon& operator*=(const Polygon& other);
	Polygon& operator^=(const T e);

	Polygon operator+(const Polygon& other) const;
	Polygon operator*(const Polygon& other) const;
	Polygon operator^(const T e) const;
		
	const VectorList& getVertices() const;
	const size_t getNumVertices() const;

	void getCCWIndices(std::vector<size_t>& indices) const;
	
	VectorList vertices;

protected:

	void topHull(std::vector<size_t>& hull) const;
	void bottomHull(std::vector<size_t>& hull) const;
	
	static bool isRightTurn2D(const Vector<T>& v1,
							  const Vector<T>& v2,
							  const Vector<T>& v3);
		
	void combineDuplicatePoints();
	void removeRedundantPoints();
	void removeRedundantPoints(const std::vector<bool>& redundant);
};

template<typename T>
std::ostream& operator<<(std::ostream& stream, const Polygon<T>& p);
	
template<typename T>
Polygon<T>::Polygon()
	: vertices() {
}

template<typename T>
Polygon<T>::Polygon(const Vector<T>& v)
	: vertices(1, v) {
}

template<typename T>
const typename Polygon<T>::VectorList& Polygon<T>::getVertices() const {
	return vertices;
}

template<typename T>
const size_t Polygon<T>::getNumVertices() const {
	return getVertices().size();
}

template<typename T>
inline Polygon<T> Polygon<T>::operator+(const Polygon<T>& other) const {
	return Polygon<T>(*this) += other;
}

template<typename T>
inline Polygon<T> Polygon<T>::operator*(const Polygon<T>& other) const {
	return Polygon<T>(*this) *= other;
}

template<typename T>
inline Polygon<T> Polygon<T>::operator^(const T e) const {
	return Polygon<T>(*this) ^= e;
}
	
template<typename T>
Polygon<T>& Polygon<T>::operator+=(const Polygon& p) {
	if (vertices.empty()) {
		*this = p;
		return *this;
	} else if (p.vertices.empty()) {
		return *this;
	} else if (p.vertices.size() == 1) {
		Vector<T> v(p.vertices.front());

		typename VectorList::iterator vertexPos =
			std::lower_bound(vertices.begin(), vertices.end(), v);

		if (vertexPos != vertices.end() and *vertexPos == v) {
			return *this;
		} else {
			vertices.insert(vertexPos, v);
		}
	} else {
		typedef typename VectorList::iterator VertexIterator;

		VectorList new_vertices(vertices.size() + p.vertices.size(),
								Vector<T>(2));
		std::merge(vertices.begin(), vertices.end(),
				   p.vertices.begin(), p.vertices.end(),
				   new_vertices.begin());

		vertices.swap(new_vertices);
		combineDuplicatePoints();
	}

	removeRedundantPoints();
	return *this;
}

template<typename T>
Polygon<T>& Polygon<T>::operator*=(const Polygon& p) {
	if (vertices.empty() or p.vertices.empty()) {
		vertices.clear();
		return *this;
	} else if (p.vertices.size() == 1) {
		Vector<T> v(p.vertices.front());
		for (size_t i = 0; i < vertices.size(); ++i) {
			vertices[i] += v;
		}
		return *this;
	} else {
		VectorList newVertices;

		newVertices.reserve(vertices.size() * p.vertices.size());

		for (size_t i = 0; i < vertices.size(); ++i) {
			for (size_t j = 0; j < p.vertices.size(); ++j) {
				newVertices.push_back(vertices[i] + p.vertices[j]);
			}
		}

		vertices.swap(newVertices);
			
		std::sort(vertices.begin(), vertices.end());
		combineDuplicatePoints();
		removeRedundantPoints();
		return *this;
	}
}
	
template<typename T>
Polygon<T>& Polygon<T>::operator^=(const T e) {
	typedef typename VectorList::iterator VectorListIter;
	for (VectorListIter v = vertices.begin(); v != vertices.end(); ++v) {
		(*v) *= e;
	}
	if (not vertices.empty() and e == 0) {
		vertices.resize(1);
	}
	return *this;
}

template<typename T>
void 
Polygon<T>::combineDuplicatePoints() {
	vertices.erase(std::unique(vertices.begin(), vertices.end()),
				   vertices.end());
}

template<typename T>
bool
Polygon<T>::isRightTurn2D(const Vector<T>& v1,
						  const Vector<T>& v2,
						  const Vector<T>& v3) {
	return ((v1[0] - v2[0]) * (v3[1] - v2[1]) -
			(v3[0] - v2[0]) * (v1[1] - v2[1])) > 0;
}

template<typename T>
void
Polygon<T>::getCCWIndices(std::vector<size_t>& indices) const {
	bottomHull(indices);
	indices.pop_back();
	topHull(indices);
	std::reverse(indices.begin(), indices.end());
	indices.pop_back();
}

template<typename T>
void
Polygon<T>::topHull(std::vector<size_t>& hull) const {
	int n = vertices.size();
	hull.push_back(0);
	hull.push_back(1);
	for (int i = 2; i < n; ++i) {
		hull.push_back(i);
		while (hull.size() > 2 and
			   not isRightTurn2D(vertices[hull[hull.size() - 3]],
								 vertices[hull[hull.size() - 2]],
								 vertices[hull[hull.size() - 1]])) {
			hull.erase(hull.end() - 2);
		}
	}
}

template<typename T>
void
Polygon<T>::bottomHull(std::vector<size_t>& hull) const {
	int n = vertices.size();
	hull.push_back(n - 1);
	hull.push_back(n - 2);
	for (int i = n - 3; i >= 0; --i) {
		hull.push_back(i);
		while (hull.size() > 2 and
			   not isRightTurn2D(vertices[hull[hull.size() - 3]],
								 vertices[hull[hull.size() - 2]],
								 vertices[hull[hull.size() - 1]])) {
			hull.erase(hull.end() - 2);
		}
	}
}

template<typename T>
void 
Polygon<T>::removeRedundantPoints() {
	int n = vertices.size();
		
	if (n <= 2) { return; }

	std::vector<bool> redundant(n, true);

	std::vector<size_t> side;
	side.reserve(n);

	topHull(side);
	for (size_t i = 0; i < side.size(); ++i) {
		redundant[side[i]] = false;
	}
	
	side.clear();

	bottomHull(side);
	for (size_t i = 0; i < side.size(); ++i) {
		redundant[side[i]] = false;
	}

	Polygon<T>::removeRedundantPoints(redundant);
}
	
template<typename T>
void 
Polygon<T>::removeRedundantPoints(const std::vector<bool>& redundant) {
	typedef typename VectorList::iterator VertexIterator;
	typedef std::vector<bool>::const_iterator RedundantIterator;
		
	VertexIterator insertV = vertices.begin(), currV = vertices.begin();
	RedundantIterator currR = redundant.begin();

	while (currV != vertices.end()) {
		if (not *currR) {
			*insertV = *currV;
			++insertV;
		}
		++currV;
		++currR;
	}
		
	vertices.erase(insertV, vertices.end());
}
		
template<typename T>
std::ostream& operator<<(std::ostream& stream, const Polygon<T>& p) {
	typedef typename Polygon<T>::VectorList::const_iterator VertexIterator;
	for (VertexIterator v = p.vertices.begin(); v != p.vertices.end(); ++v) {
		stream << *v << '\n';
	}
	return stream;
}

template<typename T>
std::istream& operator>>(std::istream& stream, Polygon<T>& p) {
	p.vertices.clear();
	Vector<T> v(2);
	while (stream >> v) {
		p.vertices.push_back(v);
	}
	return stream;
}

#endif // __POLYGON_HH__
