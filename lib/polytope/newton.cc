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

#include <Vector.h>
#include <Bitset.h>
#include <Array.h>
#include <IncidenceMatrix.h>
#include <Set.h>
#include <lrs_interface.h>
#include <cdd_interface.h>
using namespace polymake::polytope;

// #include <iostream>
// #include <ostream_wrapper.h>
// using namespace polymake;

// void printRow(const std::vector<Polytope::Newton>& v) {
// 	for (size_t j = 0; j < v.size(); ++j) {
// 		std::cerr << j << '\n';
// 		wrap(std::cerr) << v[j].getVertices() << '\n';
// 	}
// }



#include "polytope/newton.hh"

namespace polymake { namespace polytope {
	void reduceVertices(ListMatrix< Vector<Rational> >& vertices) {
		typedef ListMatrix< Vector<Rational> > VerticesType;
		std::cerr << "Creating solver...\n";
		cdd_interface::solver<Rational> s;
		std::cerr << "Finding vertices among points ("
				  << vertices.rows() << ")...\n";
		Bitset keep = s.find_vertices_among_points(vertices).first;

		std::cerr << "Removing non-vertex points...\n";
		Rows<VerticesType>& verticesRows = rows(vertices);
		Rows<VerticesType>::iterator rowIt = rows(vertices).begin();
		int rowNum = 0;
		while (rowIt != verticesRows.end()) {
			if (!keep.contains(rowNum)) {
				vertices.delete_row(rowIt++);
			} else {
				++rowIt;
			}
			++rowNum;
		}
		std::cerr << "Reducing vertices completed\n";
	}


	template <typename E, typename Matrix1, typename Matrix2> inline
	IncidenceMatrix<>
	incidence_matrix(const GenericMatrix<Matrix1,E>& R, const GenericMatrix<Matrix2,E>& C)
	{
		return IncidenceMatrix<> (R.rows(), C.rows(),
								  attach_operation(product(rows(R), rows(C), operations::mul()),
												   operations::composed11<conv<E,bool>, operations::logical_not>()).begin());
	}
		
} }

namespace polytope {

	const int Newton::POINT_FACTOR = 3;

	
	void Newton::reduce() const {
		if (reduced) {
			return;
		}

		polymake::polytope::reduceVertices(vertices);

		reduced = true;
		lastReducedSize = vertices.rows();
	}

	void Newton::reduceIfBig() const {
		if (vertices.rows() > POINT_FACTOR * lastReducedSize) {
			reduce();
		}
	}
	
	void Newton::computeFacets() const {
		if (facetsComputed) {
			return;
		}
		
		lrs_interface::solver s;
		
		typedef std::pair<Matrix<Rational>, Matrix<Rational> >  matrix_pair;
		matrix_pair mp = s.enumerate_facets(getVertices());

		facets = mp.first;
		
		facetsComputed = true;
	}

	const Newton::VerticesType& Newton::getVertices() const {
		reduce();
		return vertices;
	}

	const Newton::FacetsType& Newton::getFacets() const {
		computeFacets();
		return facets;
	}
	
	const size_t Newton::getNumVertices() const {
		return getVertices().rows();
	}
	
	const size_t Newton::getNumFacets() const {
		return getFacets().rows();
	}

	Newton& Newton::operator+=(const Newton& other) {
		vertices /= other.vertices;
		reduced = facetsComputed = false;
		return *this;
	}

	Newton& Newton::operator+=(Rational r) {
		return *this += Constant(r, vertices.cols() - 1);
	}
	
	Newton& Newton::operator*=(const Newton& other) {
		reduceIfBig();
		other.reduceIfBig();
		VerticesType newVertices(vertices.rows() * other.vertices.rows(),
								 vertices.cols(),
								 product(rows(vertices),
										 rows(other.vertices),
										 operations::add()).begin());
		cols(newVertices).front().fill(1);
		swap(vertices, newVertices);
		if (vertices.rows() > 1 && other.vertices.rows() > 1) {
			reduced = facetsComputed = false;
		}
		return *this;
	}

	Newton& Newton::operator*=(Rational r) {
		return *this *= Constant(r, vertices.cols() - 1);
	}
	
	Newton& Newton::operator^=(Rational r) {
		vertices *= r;
		cols(vertices).front().fill(1);
		return *this;
	}

	Variable::Variable(int varIndex, int numVars)
		: Newton(numVars) {
		Vector<Rational> v(unit_vector<Rational>(numVars + 1, varIndex + 1));
		v[0] = 1;
		vertices /= v;
		reduced = true;
		lastReducedSize = 1;
	}

	Constant::Constant(Rational r, int numVars)
		: Newton(numVars) {
		if (r != 0) {
			vertices /= unit_vector<Rational>(numVars + 1, 0);
			lastReducedSize = 1;
		}
		reduced = true;
	}

	// Returns an array of polytopes representing the cones of the
	// normal fan of this polytope
	Array<Poly> Newton::getNormalFan() {
		Poly np(0, ios::out | ios::in | ios::trunc, "RationalPolytope");
		np.take("VERTICES") << getVertices();

		Matrix<Rational> facets;
		np.give("FACETS") >> facets;
		
		Array< Set<int> > facetsThruVertices;
		np.give("FACETS_THRU_VERTICES") >> facetsThruVertices;

		Array<Poly> normalFan(facetsThruVertices.size());

		// Loop over the indices of the vertices of the polytope
 		for (int i = 0; i < facetsThruVertices.size(); ++i) {
			normalFan[i].init(0, ios::in | ios::out | ios::trunc,
							  "RationalPolytope");
			
 			ListMatrix< Vector<Rational> > rays;

			// For each facet thru this vertex, add the facet normal
 			typedef Entire< Set<int> >::const_iterator EntireIt;
 			for (EntireIt f = entire(facetsThruVertices[i]); !f.at_end(); ++f) {
 				Vector<Rational> facet = facets[*f];

 				// Convert to facet normal ray
 				facet *= -1;
 				facet[0] = 0;

 				// Add to rays
 				rays /= facet;
 			}
			
 			normalFan[i].take("VERTICES") << rays;
		}
		return normalFan;
	}
	
 	Array< Set<int> > Newton::getFacetsThruVertices() {
		const IncidenceMatrix<> VIF = incidence_matrix(getFacets(),
													   getVertices());
		Array< Set<int> > incidence(vertices.rows());
		for (int i = 0; i < vertices.rows(); ++i) {
			for (int j = 0; j < facets.rows(); ++j) {
				if (VIF(i, j)) {
					incidence[i].push_back(j);
				}
			}
		}

		return incidence;
 	}

	// Returns an array of polytopes representing the cones of the
	// normal fan of this polytope
	Array< Vector<Integer> > Newton::getFacetNormals() {
		computeFacets();
		
		Array<Vector<Integer> > facetNormals(facets.rows());
		
		for (int i = 0; i < facets.rows(); ++i) {
			facetNormals[i] = Vector<Integer>(facets[i].size());
			facetNormals[i][0] = 0;
			Integer scale = -1;
			for (int j = 1; j < facets[i].size(); ++j) {
				scale *= abs(denominator(facets[i][j]));
			}
			for (int j = 1; j < facets[i].size(); ++j) {
				facetNormals[i][j] = numerator(scale * facets[i][j]);
			}
			scale = gcd(facetNormals[i][1], facetNormals[i][2]);
			for (int j = 3; j < facetNormals[i].size(); ++j) {
				scale = gcd(scale, facetNormals[i][j]);
			}
			facetNormals[i] /= abs(scale);
		}

		return facetNormals;
	}

};
