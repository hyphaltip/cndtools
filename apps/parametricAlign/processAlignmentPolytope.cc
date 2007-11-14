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

Set<int> getFacetsThruPoint(const Vector<Rational>& point,
							const Matrix<Rational>& facets) {
	std::vector<int> indices;
	for (int i = 0; i < facets.rows(); ++i) {
		if (dot_product(point, facets.row(i)) == 0) {
			indices.push_back(i);
		}
	}
	return Set<int>(indices.begin(), indices.end());
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

Matrix<Rational> getFacets(const Matrix<Rational>& facets,
						   const Set<int>& indices) {
	ListMatrix< Vector<Rational> > f;
	for (Set<int>::const_iterator it = indices.begin();
		 it != indices.end(); ++it) {
		f /= facets.row(*it);
	}
	return f;
}

Matrix<Rational> computeFacets(const Matrix<Rational>& points) {
	polymake::polytope::cdd_interface::solver<Rational> solver;
	return solver.enumerate_facets(points).first;
}

bool isFeasible(const Matrix<Rational>& inequalities) {
	polymake::polytope::cdd_interface::solver<Rational> solver;
	Matrix<Rational> equations;
	try {
		Matrix<Rational> vertices = solver.enumerate_vertices(inequalities,
															  equations);
		return vertices.rows() >= inequalities.cols() - 1;
	} catch (const polymake::not_pointed& e) {
		return true;
	}
}

void makeFacets(AnnotatedPolytope& polytope) {
	if (polytope.hasSection("FACETS")) { return; }

	MatrixSection<Rational>* pointsSection;
	if (polytope.hasSection("VERTICES")) {
		polytope.getSection("VERTICES", pointsSection);
	} else if (polytope.hasSection("POINTS")) {
		polytope.getSection("POINTS", pointsSection);
	} else {
		throw std::runtime_error("Polytope has neither POINTS nor VERTICES");
	}

	Matrix<Rational> points = pointsSection->matrix;
		
	Matrix<Rational> facets = computeFacets(points);

	MatrixSection<Rational>* facetsSection = new MatrixSection<Rational>();
	facetsSection->title = "FACETS";
	facetsSection->matrix = facets;
	polytope.setSection(facetsSection);
}

void processPointsAndPointsInFacets(AnnotatedPolytope& polytope) {
	MatrixSection<Rational>* pointsSection;
	polytope.getSection("POINTS", pointsSection);
	Matrix<Rational> points = pointsSection->matrix;
	
	MatrixSection<Rational>* facetsSection;
	polytope.getSection("FACETS", facetsSection);
	Matrix<Rational> facets = facetsSection->matrix;

	SetListSection* pointsInFacetsSection;
	polytope.getSection("POINTS_IN_FACETS", pointsInFacetsSection);
	std::vector< Set<int> > pointsInFacets = pointsInFacetsSection->sets;

	// Make facets thru vertices
	std::vector< Set<int> > facetsThruPoints(points.rows());
	for (size_t i = 0; i < pointsInFacets.size(); ++i) {
		const Set<int>& pointIndices = pointsInFacets[i];
		for (Set<int>::const_iterator it = pointIndices.begin();
			 it != pointIndices.end(); ++it) {
			facetsThruPoints.at(*it).insert(static_cast<int>(i));
		}
	}
	
	int dim = points.cols() - 1;
	
	ListMatrix< Vector<Rational> > vertices;
	std::vector< Set<int> > facetsThruVertices;
	for (int i = 0; i < points.rows(); ++i) {
		Vector<Rational> point = points.row(i);
		Set<int>& indices = facetsThruPoints.at(i);
		Matrix<Rational> facetInequalities = getFacets(facets, indices);
		if (rank(facetInequalities) == dim) {
			vertices /= point;
			facetsThruVertices.push_back(indices);
		}
	}

	MatrixSection<Rational>* vertsSection = new MatrixSection<Rational>();
	vertsSection->title = "VERTICES";
	vertsSection->matrix = vertices;
	polytope.setSection(vertsSection);
		
	SetListSection* facetsThruVerticesSection = new SetListSection();
	facetsThruVerticesSection->title = "FACETS_THRU_VERTICES";
	facetsThruVerticesSection->sets = facetsThruVertices;
	polytope.setSection(facetsThruVerticesSection);
}

void makeVertices(AnnotatedPolytope& polytope) {
	if (polytope.hasSection("VERTICES") and
		polytope.hasSection("FACETS_THRU_VERTICES")) {
		return;
	}

	if (polytope.hasSection("POINTS") and
		polytope.hasSection("POINTS_IN_FACETS")) {
		processPointsAndPointsInFacets(polytope);
		return;
	}
	
	MatrixSection<Rational>* pointsSection;
	bool computeVertices = true;
	if (polytope.hasSection("VERTICES")) {
		polytope.getSection("VERTICES", pointsSection);
		computeVertices = false;
	} else if (polytope.hasSection("POINTS")) {
		polytope.getSection("POINTS", pointsSection);
	} else {
		throw std::runtime_error("Polytope has neither POINTS nor VERTICES");
	}

	Matrix<Rational> points = pointsSection->matrix;
	
	MatrixSection<Rational>* facetsSection;
	polytope.getSection("FACETS", facetsSection);
	Matrix<Rational> facets = facetsSection->matrix;
		
	int dim = points.cols() - 1;
	
	ListMatrix< Vector<Rational> > vertices;
	std::vector< Set<int> > facetsThruVertices;
	for (int i = 0; i < points.rows(); ++i) {
		Vector<Rational> point = points.row(i);
		Set<int> indices = getFacetsThruPoint(point, facets);
		Matrix<Rational> facetInequalities = getFacets(facets, indices);
		if (computeVertices) {
			if (rank(facetInequalities) == dim) {
				vertices /= point;
				facetsThruVertices.push_back(indices);
			}
		} else {
			facetsThruVertices.push_back(indices);
		}			
	}

	if (computeVertices) {
		MatrixSection<Rational>* vertsSection = new MatrixSection<Rational>();
		vertsSection->title = "VERTICES";
		vertsSection->matrix = vertices;
		polytope.setSection(vertsSection);
	}
		
	SetListSection* facetsThruVerticesSection = new SetListSection();
	facetsThruVerticesSection->title = "FACETS_THRU_VERTICES";
	facetsThruVerticesSection->sets = facetsThruVertices;
	polytope.setSection(facetsThruVerticesSection);
}

void makeFacetNormals(AnnotatedPolytope& polytope) {
	if (not polytope.hasSection("FACETS")) {
		makeFacets(polytope);
	}

	MatrixSection<Rational>* facetsSection;
	polytope.getSection("FACETS", facetsSection);
	Matrix<Rational> facets = facetsSection->matrix;

	Matrix<Integer> facetNormals(facets.rows(), facets.cols());
	for (int i = 0; i < facetNormals.rows(); ++i) {
		Vector<Rational> facet = facets.row(i);
		facet[0] = 0;
		facetNormals.row(i) = integerize(facet);
		facetNormals.row(i) *= -1;
	}

	MatrixSection<Integer>* facetNormalsSection = new MatrixSection<Integer>();
	facetNormalsSection->title = "FACET_NORMALS";
	facetNormalsSection->matrix = facetNormals;
	polytope.setSection(facetNormalsSection);
}

void makeReasonableVertices(AnnotatedPolytope& polytope,
						    Matrix<Rational>& inequalities) {
	MatrixSection<Integer>* facetNormalsSection;
	polytope.getSection("FACET_NORMALS", facetNormalsSection);
	Matrix<Integer> facetNormals = facetNormalsSection->matrix;

	SetListSection* facetsThruVerticesSection;
	polytope.getSection("FACETS_THRU_VERTICES", facetsThruVerticesSection);
	std::vector< Set<int> > facetsThruVertices = 
		dynamic_cast<SetListSection*>(facetsThruVerticesSection)->sets;

	Matrix<Rational> isReasonable(facetsThruVertices.size(), 1);
	for (size_t i = 0; i < facetsThruVertices.size(); ++i) {
		Matrix<Integer> rays = getFacetNormals(facetNormals,
											   facetsThruVertices[i]);
		Matrix<Rational> ratRays = convert_to<Rational>(rays);
		Matrix<Rational> cone = computeFacets(ratRays);
		cone /= inequalities;

		isReasonable(i, 0) = isFeasible(cone);
	}
	
	MatrixSection<Rational>* reasonableVerticesSection = new MatrixSection<Rational>();
	reasonableVerticesSection->title = "REASONABLE_VERTICES";
	reasonableVerticesSection->matrix = isReasonable;
	polytope.setSection(reasonableVerticesSection);
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	std::string reasonableFilename;

	// Set up option parser
	util::options::Parser parser("", "");
	parser.addStoreOpt('r', "reasonable",
					   "Polymake file containing inequalities for reasonable "
					   "parameter values", reasonableFilename);
	parser.parse(argv, argv + argc);

	try {
		Matrix<Rational> inequalities;
		if (not reasonableFilename.empty()) {
			InputFileStream reasonableFile(reasonableFilename);
			AnnotatedPolytope reasonablePoly;
			reasonableFile >> reasonablePoly;
			if (not reasonablePoly.hasSection("INEQUALITIES")) {
				throw std::runtime_error("Reasonable file does not have "
										 "section INEQUALITIES");
			}
			MatrixSection<Rational>* inequalitiesSection;
			reasonablePoly.getSection("INEQUALITIES", inequalitiesSection);
			inequalities = inequalitiesSection->matrix;
		}
		
		AnnotatedPolytope polytope;
		std::cin >> polytope;

		makeFacets(polytope);
		makeVertices(polytope);
		makeFacetNormals(polytope);

		if (inequalities.rows() != 0) {
			makeReasonableVertices(polytope, inequalities);
		}
		
		std::cout << polytope;
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
