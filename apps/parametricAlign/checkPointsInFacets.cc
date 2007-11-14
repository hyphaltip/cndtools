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
using pm::incl;
using pm::rank;

Rational dot_product(const Vector<Rational>& v1, const Vector<Rational>& v2) {
	Rational r = 0;
	for (int i = 0; i < v1.size(); ++i) {
		r += v1[i] * v2[i];
	}
	return r;
}

Set<int> getPointsInFacet(const Vector<Rational>& facet,
						  const Matrix<Rational>& points) {
	std::vector<int> indices;
	for (int i = 0; i < points.rows(); ++i) {
		if (dot_product(facet, points.row(i)) == 0) {
			indices.push_back(i);
		}
	}
	return Set<int>(indices.begin(), indices.end());
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

		MatrixSection<Rational>* pointsSection;
		polytope.getSection("POINTS", pointsSection);
		Matrix<Rational> points = pointsSection->matrix;

		MatrixSection<Rational>* facetsSection;
		polytope.getSection("FACETS", facetsSection);
		Matrix<Rational> facets = facetsSection->matrix;

		SetListSection* pointsInFacetsSection;
		polytope.getSection("POINTS_IN_FACETS", pointsInFacetsSection);
		std::vector< Set<int> >pointsInFacets = pointsInFacetsSection->sets;

		for (int i = 0; i < facets.rows(); ++i) {
			Set<int> s = getPointsInFacet(facets.row(i), points);
			if (s != pointsInFacets[i]) {
				std::cerr << "Points in facets error for facet " << i << '\n';
				wrap(std::cerr) << s << '\n';
				wrap(std::cerr) << pointsInFacets[i] << '\n';
				
			}
		}
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
