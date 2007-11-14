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

#include "polytope/Cone.hh"

namespace polytope {

	polymake::polytope::cdd_interface::solver<Rational> Cone::CDDSolver;

	Cone::Ray Cone::ratToInt(const Vector<Rational>& ratRay) {
		Ray intRay(ratRay.size());

		// Clear denominators
		Integer scale = 1;
		for (size_t i = 1; i < ratRay.size(); ++i) {
			scale *= abs(denominator(ratRay[i]));
		}
		for (size_t i = 1; i < ratRay.size(); ++i) {
			intRay[i] = numerator(scale * ratRay[i]);
		}

		// Reduce by greatest common divisor
		scale = gcd(intRay[1], intRay[2]);
		for (size_t i = 3; i < ratRay.size(); ++i) {
			scale = gcd(scale, intRay[i]);
		}
		scale = abs(scale);
		intRay /= scale;

		return intRay;
	}
	
	Cone::Cone(const RayList& rays)
		: rays(rays),
		  facets(),
		  reduced(false) {
		computeFacets();
		this->rays.clear();
	}
	
	void Cone::addInequality(const Inequality& i) {
		facets.push_back(i);
		reduced = false;
	}
				
	size_t Cone::getNumFacets() {
		reduce();
		return facets.size();
	}
	
	Cone::Facet Cone::getFacet(size_t i) {
		reduce();
		return facets.at(i);
	}
	
	size_t Cone::getNumRays() {
		reduce();
		return rays.size();
	}
	
	Cone::Ray Cone::getRay(size_t i) {
		reduce();
		return rays.at(i);
	}

	void Cone::reduce() {
		if (reduced) { return; }
		computeRays();
		computeFacets();
		reduced = true;
	}
		
	void Cone::computeRays() {
		rays.clear();
		
		if (facets.empty()) { return; }

		size_t dim = facets.front().size() - 1;
		polymake::Matrix<Rational> inequalities(facets.size(), dim + 1);
		polymake::Matrix<Rational> equations; //(0, dim + 1);

		for (size_t i = 0; i < facets.size(); ++i) {
			for (size_t j = 0; j < dim + 1; ++j) {
				inequalities(i, j) = facets.at(i)[j];
			}
		}

		polymake::Matrix<Rational> ratRays =
			CDDSolver.enumerate_vertices(inequalities, equations);

		for (int i = 0; i < ratRays.rows(); ++i) {
			if (ratRays(i, 0) != 0) { continue; } // This must be the origin
			Vector<Rational> ratRay(ratRays.cols());
			for (int j = 0; j < ratRays.cols(); ++j) {
				ratRay[j] = ratRays(i, j);
			}
			rays.push_back(ratToInt(ratRay));
		}
	}
		
	void Cone::computeFacets() {
		facets.clear();
		
		if (rays.empty()) { return; }
		
		size_t dim = rays.front().size() - 1;
		polymake::Matrix<Rational> points(rays.size(), dim + 1);
										  
		for (size_t i = 0; i < rays.size(); ++i) {
			for (size_t j = 0; j < dim + 1; ++j) {
				points(i, j) = rays.at(i)[j];
			}
		}

		polymake::Matrix<Rational> inequalities =
			CDDSolver.enumerate_facets(points).first;

		for (int i = 0; i < inequalities.rows(); ++i) {
			Vector<Rational> ratFacet(inequalities.cols());
			for (int j = 0; j < inequalities.cols(); ++j) {
				ratFacet[j] = inequalities(i, j);
			}
			facets.push_back(ratToInt(ratFacet));
		}
	}

}
