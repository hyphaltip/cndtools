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

#ifndef __POLYTOPE_CONE_HH__
#define __POLYTOPE_CONE_HH__

#include <vector>

#include "polytope/Vector.hh"

// Polymake includes
#include <Rational.h>
#include <Integer.h>
#include <Matrix.h>
#include <Bitset.h>
#include <cdd_interface.h>

namespace polytope {

    class Cone {
	public:
		typedef Vector<Integer> Ray;
		typedef std::vector<Ray> RayList;
		typedef Vector<Integer> Inequality;
		typedef Inequality Facet;
		typedef std::vector<Facet> FacetList;
		
		Cone(const RayList& rays);

		void addInequality(const Inequality& i);
		
		size_t getNumFacets();
		Facet getFacet(size_t i);

		size_t getNumRays();
		Ray getRay(size_t i);
		
	private:
		void reduce();
		void computeRays();
		void reduceRays();
		void computeFacets();

		static Ray ratToInt(const Vector<Rational>& ratRay);

		static polymake::polytope::cdd_interface::solver<Rational> CDDSolver;
		
		RayList rays;
		FacetList facets;
		bool reduced;
    };

}

#endif // __POLYTOPE_CONE_HH__
