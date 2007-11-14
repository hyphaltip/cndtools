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

#ifndef __POLYTOPE_FORMATS_POLYMAKE_OUTPUTSTREAM_HH__
#define __POLYTOPE_FORMATS_POLYMAKE_OUTPUTSTREAM_HH__

#include <iosfwd>

#include "polytope/Polytope.hh"

namespace polytope { namespace formats { namespace polymake {

    class OutputStream {
	public:
		OutputStream(std::ostream& strm);

		template<typename T>
		OutputStream& operator<<(const Polytope<T>& p);

	private:
		std::ostream& strm;
    };

	inline OutputStream::OutputStream(std::ostream& strm) : strm(strm) {}

	template<typename T>
	OutputStream& OutputStream::operator<<(const Polytope<T>& p) {
		const typename Polytope<T>::VectorList& vertices = p.getVertices();
		strm << "VERTICES" << '\n';
		for (size_t i = 0; i < vertices.size(); ++i) {
			strm << "1\t" << vertices[i] << '\n';
		}
		strm << '\n';
// 		const typename Polytope<T>::CountList& counts = p.getVertexCounts();
// 		strm << "VERTEX_COUNTS" << '\n';
// 		for (size_t i = 0; i < counts.size(); ++i) {
// 			strm << counts[i] << '\n';
// 		}
// 		strm << '\n';
// 		const typename Polytope<T>::FacetList& facets = p.getFacets();
// 		strm << "FACETS" << '\n';
// 		for (size_t i = 0; i < facets.size(); ++i) {
// 			strm << facets[i] << '\n';
// 		}
// 		strm << '\n';
// 		const typename Polytope<T>::FacetNormalList facetNormals = p.getFacetNormals();
// 		strm << "FACET_NORMALS" << '\n';
// 		for (size_t i = 0; i < facetNormals.size(); ++i) {
// 			strm << facetNormals[i] << '\n';
// 		}
// 		strm << '\n';
// 		const typename Polytope<T>::FacetsThruVertexList facetsThruVertices = p.getFacetsThruVertices();
// 		strm << "FACETS_THRU_VERTICES" << '\n';
// 		for (size_t i = 0; i < facetsThruVertices.size(); ++i) {
// 			strm << '{';
// 			typedef typename Polytope<T>::FacetsThruVertex::const_iterator FacetIndexIterator;
// 			for (FacetIndexIterator fi = facetsThruVertices[i].begin();
// 				 fi != facetsThruVertices[i].end(); ++fi) {
// 				if (fi != facetsThruVertices[i].begin()) { strm << ' '; }
// 				strm << *fi;
// 			}
// 			strm << "}\n";
// 		}
// 		strm << '\n';
 		return *this;
	}
		
} } }

#endif // __POLYTOPE_FORMATS_POLYMAKE_OUTPUTSTREAM_HH__
