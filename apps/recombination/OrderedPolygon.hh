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

#ifndef __ORDEREDPOLYGON_HH__
#define __ORDEREDPOLYGON_HH__

#include <vector>
#include <iosfwd>

#include "Polygon.hh"

template<typename T>
class OrderedPolygon {
public:
	typedef typename std::vector< Vector<T> > VectorList;

	OrderedPolygon() {
	}
	
	OrderedPolygon(const Polygon<T>& p) {
		std::vector<size_t> ordered_indices;
		p.getCCWIndices(ordered_indices);
		for (size_t i = 0; i < ordered_indices.size(); ++i) {
			vertices.push_back(p.vertices[ordered_indices[i]]);
		}
		makeRays();
	}

	size_t getNumVertices() const {
		return vertices.size();
	}
	
	template<typename coord_type>
	size_t getVertexNum(const Vector<coord_type>& ray) const {
		double angle = rayToAngle(ray);
		std::vector<double>::const_iterator it =
			std::lower_bound(rays.begin(), rays.end(), angle);
		return std::distance(rays.begin(), it);
	}

	Vector<T> ccwRay(size_t i) const {
		size_t j = (i == vertices.size() - 1 ? 0 : i + 1);
		Vector<T> halfspace = vertices[i] - vertices[j];
		return makeRay(-halfspace[1], halfspace[0]);
	}

	Vector<T> cwRay(size_t i) const {
		return ccwRay(i == 0 ? vertices.size() - 1 : i - 1);
	}

	Vector<T> vertexNormal(size_t i) const {
		return ccwRay(i) + cwRay(i);
	}
	
	std::pair< Vector<T>, Vector<T> >getCone(size_t i) {
		return std::make_pair(ccwRay(i), cwRay(i));
	}

	template<typename coord_type>
	static double rayToAngle(const Vector<coord_type>& v) {
		double angle = std::atan2(static_cast<double>(v[1]),
								  static_cast<double>(v[0]));
		if (angle <= 0.0) { angle += 2 * M_PI; }
		return angle;
	}

	template<typename coord_type>
	static Vector<coord_type> makeRay(coord_type x, coord_type y) {
		Vector<coord_type> ray(2);
		ray[0] = x;
		ray[1] = y;
		return ray;
	}

	friend std::ostream& operator<<(std::ostream& stream,
									const OrderedPolygon<T>& p) {
		for (size_t i = 0; i < p.vertices.size(); ++i) {
			stream << p.vertices[i] << '\n';
		}
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream,
									OrderedPolygon<T>& p) {
		p.vertices.clear();
		Vector<T> v(2);
		while (stream >> v) {
			p.vertices.push_back(v);
		}
		p.makeRays();
		return stream;
	}
	
protected:

	void makeRays() {
		rays.clear();
		for (size_t i = 0; i < vertices.size(); ++i) {
			rays.push_back(rayToAngle(ccwRay(i)));
		}
	}
	
	VectorList vertices;
	std::vector<double> rays;
};

#endif // __ORDEREDPOLYGON_HH__
