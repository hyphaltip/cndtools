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

#include <iostream>

#include "util/options.hh"
#include "util/io/line/InputStream.hh"
#include "filesystem.hh"
using namespace filesystem;

#include "AnnotatedPolytope.hh"

#include <linalg.h>
#include <Rational.h>
#include <Matrix.h>
#include <ListMatrix.h>
#include <Vector.h>
#include <Set.h>
#include <cdd_interface.h>
#include <ostream_wrapper.h>

using polymake::ListMatrix;
using polymake::Matrix;
using polymake::Vector;
using polymake::Set;
using polymake::convert_to;
using polymake::wrap;
using pm::rank;

Rational dot_product(const Vector<Rational>& v1, const Vector<Rational>& v2) {
	Rational r = 0;
	for (int i = 0; i < v1.size(); ++i) {
		r += v1[i] * v2[i];
	}
	return r;
}

Vector<Integer> integerize(const Vector<Rational>& v) {
	Integer scale = 1;
	Vector<Integer> iv(v.size());
	for (int i = 0; i < v.size(); ++i) {
		scale *= abs(denominator(v[i]));
	}
	for (int i = 0; i < v.size(); ++i) {
		iv[i] = numerator(v[i]) * (scale / denominator(v[i]));
	}
	scale = gcd(iv[1], iv[2]);
	for (int i = 3; i < iv.size(); ++i) {
		scale = gcd(scale, iv[i]);
	}
	scale = abs(scale);
	iv /= scale;

	return iv;
}

Matrix<Integer> getFacetNormals(const Matrix<Integer>& facetNormals,
								 const Set<int>& indices) {
	ListMatrix< Vector<Integer> > f;
	for (Set<int>::const_iterator it = indices.begin();
		 it != indices.end(); ++it) {
		f /= facetNormals.row(*it);
	}
	return f;
}

Matrix<Rational> computeFacets(const Matrix<Rational>& points) {
	polymake::polytope::cdd_interface::solver<Rational> solver;
	return solver.enumerate_facets(points).first;
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments

	// Set up option parser
	util::options::Parser parser("", "");
	parser.parse(argv, argv + argc);

	try {
		AnnotatedPolytope polytope;
		std::cin >> polytope;
		
		MatrixSection<Integer>* facetNormalsSection;
		polytope.getSection("FACET_NORMALS", facetNormalsSection);
		Matrix<Integer> facetNormals = facetNormalsSection->matrix;
		
		SetListSection* facetsThruVerticesSection;
		polytope.getSection("FACETS_THRU_VERTICES", facetsThruVerticesSection);
		std::vector< Set<int> > facetsThruVertices = 
			dynamic_cast<SetListSection*>(facetsThruVerticesSection)->sets;
		
		for (size_t i = 0; i < facetsThruVertices.size(); ++i) {
			Matrix<Integer> rays = getFacetNormals(facetNormals,
												   facetsThruVertices[i]);

			Matrix<Rational> ratRays = convert_to<Rational>(rays);
			Matrix<Rational> cone = computeFacets(ratRays);
			
			for (int j = 0; j < cone.rows(); ++j) {
				Vector<Rational> ineqRat = cone.row(j);
				Vector<Integer> ineq = integerize(ineqRat);
				for (int k = 0; k < ineq.size(); ++k) {
					std::cout << ineq[k] << ' ';
				}
			}
			
			std::cout << '\n';
		}
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
