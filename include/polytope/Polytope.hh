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

#ifndef __POLYTOPE_POLYTOPE_HH__
#define __POLYTOPE_POLYTOPE_HH__

#include <vector>
#include <set>
#include <algorithm>
#include <iosfwd>

#include "polytope/Vector.hh"
#include "polytope/Cone.hh"

// Polymake includes
#include <Rational.h>
#include <Integer.h>
#include <Matrix.h>
#include <Bitset.h>
#include <lrs_interface.h>
#include <cdd_interface.h>

// Qhull includes
extern "C" {
#include "qhull/qhull_a.h"
}

namespace polytope {

	template<typename T>
    class Polytope {
	public:
		typedef Vector<Rational> Facet;
		typedef Vector<Integer> Ray;
		typedef typename std::vector< Vector<T> > VectorList;
		typedef std::vector<Facet> FacetList;
		typedef std::vector<Ray> FacetNormalList;
		typedef std::set<size_t> FacetsThruVertex;
		typedef std::vector<FacetsThruVertex> FacetsThruVertexList;
		typedef std::vector<Cone> NormalFan;

		Polytope(size_t ambientDim = 0);
		Polytope(const Vector<T>& v);
		Polytope(const VectorList& points,
				 bool reduced = false);
		
		Polytope& operator+=(const Polytope& other);
		Polytope& operator*=(const Polytope& other);
		Polytope& operator^=(const T e);

		Polytope operator+(const Polytope& other) const;
		Polytope operator*(const Polytope& other) const;
		Polytope operator^(const T e) const;
		
		size_t getAmbientDim() const;

		const VectorList& getVertices() const;
		const size_t getNumVertices() const;

		const FacetList& getFacets() const;
		const size_t getNumFacets() const;

		FacetsThruVertexList getFacetsThruVertices() const;

		Ray getFacetNormal(size_t i) const;
		FacetNormalList getFacetNormals() const;

		NormalFan getNormalFan() const;
		
		template<typename C>
		friend std::ostream& operator<<(std::ostream& strm,
										const Polytope<C>& p);

		static void setQhullLogFile(const std::string& filename);

	protected:
		static bool isRightTurn2D(const Vector<T>& v1,
								  const Vector<T>& v2,
								  const Vector<T>& v3);
		
		static bool isFacetThruVertex(const Facet& f, const Vector<T>& v);
		
		void combineDuplicatePoints() const;
		
		void removeRedundantPoints() const;
		void removeRedundantPoints2D() const;
		void removeRedundantPoints(const std::vector<bool>& redundant) const;
		void ensureReduced() const;

		bool redundantQhull(std::vector<bool>& redundant) const;
		template<typename Solver>
		bool redundantPolymake(std::vector<bool>& redundant, Solver& s) const;
		bool redundantCDD(std::vector<bool>& redundant) const;
		bool redundantLRS(std::vector<bool>& redundant) const;

		void computeFacets() const;

		template<typename Solver>
		void computeFacetsPolymake(Solver& s) const;
		void computeFacetsCDD() const;
		void computeFacetsLRS() const;

		void findValidDims(std::vector<size_t>& validDims) const;
		
		polymake::Matrix<Rational> polymakePoints() const;

		static polymake::polytope::lrs_interface::solver LRSSolver;
		static polymake::polytope::cdd_interface::solver<Rational> CDDSolver;
		static FILE* qhullLogFile;

		size_t ambientDim;
		mutable VectorList vertices;
		mutable FacetList facets;
		mutable bool reduced;
		mutable bool reducedEnsured;
		mutable bool facetsComputed;
    };

	template<typename T>
	polymake::polytope::lrs_interface::solver Polytope<T>::LRSSolver;

	template<typename T>
	polymake::polytope::cdd_interface::solver<Rational> Polytope<T>::CDDSolver;

	template<typename T>
	FILE* Polytope<T>::qhullLogFile = stderr;

	template<typename T>
	void Polytope<T>::setQhullLogFile(const std::string& filename) {
		if (qhullLogFile != stderr) { fclose(qhullLogFile); }
		qhullLogFile = fopen(filename.c_str(), "w");
	}
	
	template<typename T>
	size_t
	Polytope<T>::getAmbientDim() const {
		return ambientDim;
	}

	template<typename T>
	Polytope<T>::Polytope(size_t ambientDim)
		: ambientDim(ambientDim),
		  vertices(),
		  facets(),
		  reduced(true),
		  reducedEnsured(true),
		  facetsComputed(true) {}
	
	template<typename T>
	Polytope<T>::Polytope(const Vector<T>& v)
		: ambientDim(v.size()),
		  vertices(1, v),
		  facets(),
		  reduced(true),
		  reducedEnsured(true),
		  facetsComputed(false) {}

	template<typename T>
	Polytope<T>::Polytope(const VectorList& points, bool reduced)
		: ambientDim(points.front().size()),
		  vertices(points),
		  facets(),
		  reduced(reduced),
		  reducedEnsured(reduced),
		  facetsComputed(false) {
		  std::sort(vertices.begin(), vertices.end());
	}

	template<typename T>
	const typename Polytope<T>::VectorList& Polytope<T>::getVertices() const {
		ensureReduced();
		return vertices;
	}

	template<typename T>
	const size_t Polytope<T>::getNumVertices() const {
		return getVertices().size();
	}

	template<typename T>
	const typename Polytope<T>::FacetList& Polytope<T>::getFacets() const {
		computeFacets();
		return facets;
	}
	
	template<typename T>
	bool Polytope<T>::isFacetThruVertex(const Facet& f,
										const Vector<T>& v) {
		Rational sum = f[0];
		for (size_t i = 0; i < v.size(); ++i) {
			sum += (f[i + 1] * v[i]);
		}
		return (sum == 0);
	}

	template<typename T>
	typename Polytope<T>::FacetsThruVertexList
	Polytope<T>::getFacetsThruVertices() const {
		ensureReduced();
		computeFacets();

		FacetsThruVertexList list(vertices.size());
		
		for (size_t i = 0; i < vertices.size(); ++i) {
			for (size_t j = 0; j < facets.size(); ++j) {
				if (isFacetThruVertex(facets[j], vertices[i])) {
					list[i].insert(j);
				}
			}
		}

		return list;
	}

	template<typename T>
	typename Polytope<T>::FacetNormalList Polytope<T>::getFacetNormals() const {
		computeFacets();
		FacetNormalList facetNormals;

		for (size_t i = 0; i < facets.size(); ++i) {
			facetNormals.push_back(getFacetNormal(i));
		}
		
		return facetNormals;
	}

	template<typename T>
	typename Polytope<T>::Ray
	Polytope<T>::getFacetNormal(size_t i) const {
		computeFacets();

		Ray facetNormal(facets[i].size());
		facetNormal[0] = 0;
		Integer scale = -1;
		for (size_t j = 1; j < facets[i].size(); ++j) {
			scale *= abs(denominator(facets[i][j]));
		}
		for (size_t j = 1; j < facets[i].size(); ++j) {
			facetNormal[j] = numerator(scale * facets[i][j]);
		}
		scale = gcd(facetNormal[1], facetNormal[2]);
		for (size_t j = 3; j < facets[i].size(); ++j) {
			scale = gcd(scale, facetNormal[j]);
		}
		scale = abs(scale);
		facetNormal /= scale;

		return facetNormal;
	}

	template<typename T>
	typename Polytope<T>::NormalFan
	Polytope<T>::getNormalFan() const {
		computeFacets();

		NormalFan normalFan;
		FacetsThruVertexList facetsThruVertices = getFacetsThruVertices();
		FacetNormalList facetNormals = getFacetNormals();

		for (size_t i = 0; i < getNumVertices(); ++i) {
			FacetNormalList rays;
			typedef FacetsThruVertex::const_iterator FacetIndexIterator;
			for (FacetIndexIterator fi = facetsThruVertices[i].begin();
				 fi != facetsThruVertices[i].end(); ++fi) {
				rays.push_back(facetNormals.at(*fi));
			}
			normalFan.push_back(Cone(rays));
		}

		return normalFan;
	}
	
	template<typename T>
	void Polytope<T>::computeFacets() const {
		if (facetsComputed) { return; }

		removeRedundantPoints();

		computeFacetsCDD();

		facetsComputed = true;
	}

	template<typename T>
	template<typename Solver>
	void Polytope<T>::computeFacetsPolymake(Solver& s) const {
		polymake::Matrix<Rational> points = polymakePoints();
		polymake::Matrix<Rational> f = s.enumerate_facets(points).first;

		facets.clear();
		facets.resize(f.rows(), Vector<Rational>(f.cols()));
		for (int i = 0; i < f.rows(); ++i) {
			for (int j = 0; j < f.cols(); ++j) {
				facets[i][j] = f(i, j);
			}
		}
	}

	template<typename T>
	void Polytope<T>::computeFacetsCDD() const {
		computeFacetsPolymake(CDDSolver);
	}

	template<typename T>
	void Polytope<T>::computeFacetsLRS() const {
		computeFacetsPolymake(LRSSolver);
	}

	template<typename T>
	const size_t Polytope<T>::getNumFacets() const {
		return getFacets().size();
	}
	
	template<typename T>
	inline Polytope<T> Polytope<T>::operator+(const Polytope<T>& other) const {
		return Polytope<T>(*this) += other;
	}

	template<typename T>
	inline Polytope<T> Polytope<T>::operator*(const Polytope<T>& other) const {
		return Polytope<T>(*this) *= other;
	}

	template<typename T>
	inline Polytope<T> Polytope<T>::operator^(const T e) const {
		return Polytope<T>(*this) ^= e;
	}
	
	template<typename T>
	Polytope<T>& Polytope<T>::operator+=(const Polytope& p) {
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

			VectorList newVertices;

			newVertices.reserve(vertices.size() + p.vertices.size());

			VertexIterator v1 = vertices.begin(), v2 = p.vertices.begin();
			while (v1 != vertices.end() or v2 != p.vertices.end()) {
				if (v2 == p.vertices.end() or
					(v1 != vertices.end() and (*v1) < (*v2))) {
					newVertices.push_back(*v1);
					++v1;
				} else{
					newVertices.push_back(*v2);
					++v2;
				}
			}
			
			vertices.swap(newVertices);
			combineDuplicatePoints();
		}

		reducedEnsured = reduced = false;
		removeRedundantPoints();
		return *this;
	}

	template<typename T>
	Polytope<T>& Polytope<T>::operator*=(const Polytope& p) {
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

			reducedEnsured = reduced = false;
			removeRedundantPoints();
			return *this;
		}
	}
	
	template<typename T>
	Polytope<T>& Polytope<T>::operator^=(const T e) {
		for (typename VectorList::iterator v = vertices.begin();
			 v != vertices.end(); ++v) {
			(*v) *= e;
		}
		if (not vertices.empty() and e == 0) {
			vertices.resize(1);
		}
		return *this;
	}

	template<typename T>
	void 
	Polytope<T>::combineDuplicatePoints() const {
		if (vertices.empty()) { return; }
		
		typedef typename VectorList::iterator VertexIterator;

		VertexIterator v = vertices.begin() + 1, lastV = vertices.begin();

		while (v != vertices.end()) {
			if (not (*v == *lastV)) {
				++lastV;
				(*lastV) = (*v);
			}
			++v;
		}

		vertices.erase(lastV + 1, vertices.end());
	}

	template<typename T>
	void
	Polytope<T>::ensureReduced() const {
		if (reducedEnsured) { return; }

		if (vertices.size() > 2) {
			std::vector<bool> redundant(vertices.size(), false);
			redundantCDD(redundant);
			removeRedundantPoints(redundant);
		}

		reducedEnsured = reduced = true;
	}
	
	template<typename T>
	void 
	Polytope<T>::removeRedundantPoints() const {
		if (getAmbientDim() == 2) {
			removeRedundantPoints2D();
			return;
		}
		if (reduced) { return; }
		if (vertices.size() < 2) { reduced = true; return; }

// 		std::cerr << "Finding vertices among points ("
// 				  << vertices.size() << ")...\n";
		
		std::vector<bool> redundant(vertices.size(), false);
		
		bool success = (vertices.size() > ambientDim and
						redundantQhull(redundant));
		if (not success) {
			redundantCDD(redundant);
			reducedEnsured = true;
		}

		removeRedundantPoints(redundant);
		
		reduced = true;
	}

	template<typename T>
	bool
	Polytope<T>::isRightTurn2D(const Vector<T>& v1,
							   const Vector<T>& v2,
							   const Vector<T>& v3) {
		return ((v1[0] - v2[0]) * (v3[1] - v2[1]) -
				(v3[0] - v2[0]) * (v1[1] - v2[1])) > 0;
	}
	
	template<typename T>
	void 
	Polytope<T>::removeRedundantPoints2D() const {
		int n = vertices.size();
		
		if (reduced) { return; }
		if (n <= 2) { reduced = true; return; }

		std::vector<bool> redundant(n, true);

		std::vector<int> side;
		side.reserve(n);

		// Compute side hull
		side.push_back(0);
		side.push_back(1);
		for (int i = 2; i < n; ++i) {
			side.push_back(i);
			while (side.size() > 2 and
				   not isRightTurn2D(vertices[side[side.size() - 3]],
									 vertices[side[side.size() - 2]],
									 vertices[side[side.size() - 1]])) {
				side.erase(side.end() - 2);
			}
		}

		for (size_t i = 0; i < side.size(); ++i) {
			redundant[side[i]] = false;
		}

		side.clear();
		
		// Compute side hull
		side.push_back(n - 1);
		side.push_back(n - 2);
		for (int i = n - 3; i >= 0; --i) {
			side.push_back(i);
			while (side.size() > 2 and
				   not isRightTurn2D(vertices[side[side.size() - 3]],
									 vertices[side[side.size() - 2]],
									 vertices[side[side.size() - 1]])) {
				side.erase(side.end() - 2);
			}
		}

		for (size_t i = 0; i < side.size(); ++i) {
			redundant[side[i]] = false;
		}

		Polytope<T>::removeRedundantPoints(redundant);
		reducedEnsured = true;
		reduced = true;
	}
	
	template<typename T>
	void 
	Polytope<T>::removeRedundantPoints(const std::vector<bool>& redundant) const {
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
	void
	Polytope<T>::findValidDims(std::vector<size_t>& validDims) const {
		if (vertices.empty()) { return; }
		for (size_t d = 0; d < ambientDim; ++d) {
			int coord = vertices.front()[d];
			for (size_t i = 0; i < vertices.size(); ++i) {
				if (vertices[i][d] != coord) {
					validDims.push_back(d);
					break;
				}
			}
		}
	}
	
	template<typename T>
	bool Polytope<T>::redundantQhull(std::vector<bool>& redundant) const {
		int dim;                  /* dimension of points */
		int numpoints;            /* number of points */
		boolT ismalloc;           /* True if qhull should free points in qh_freeqhull() or reallocation */ 
		char flags[] = "qhull ";  /* option flags for qhull, see qh_opt.htm */
		FILE *outfile = NULL;     /* output from qh_produce_output()                    
								     use NULL to skip qh_produce_output() */ 
		FILE *errfile = qhullLogFile;   /* error messages from qhull code */ 
		int exitcode;             /* 0 if no error from qhull */
		int curlong, totlong;     /* memory remaining after qh_memfreeshort */
		vertexT *vertex;
		
		/* initialize dim, numpoints, points[], ismalloc here */
		std::vector<size_t> validDims;
		findValidDims(validDims);
		ismalloc = false;
		dim = validDims.size();
		numpoints = vertices.size();
		std::vector<coordT> points(numpoints * dim);
		for (int i = 0; i < numpoints; ++i) {
			for (size_t j = 0; j < validDims.size(); ++j) {
				points[i * dim + j] = vertices[i][validDims[j]];
			}
		}
		
		exitcode = qh_new_qhull (dim, numpoints, &points[0], ismalloc,
								 flags, outfile, errfile);

		if (!exitcode) { /* if no error */ 
			std::fill(redundant.begin(), redundant.end(), true);
			FORALLvertices {
				size_t pointIndex = (vertex->point - &points[0]) / dim;
				redundant.at(pointIndex) = false;
			}
		}
		qh_freeqhull(not qh_ALL);  
		qh_memfreeshort(&curlong, &totlong);
		if (curlong or totlong)
			fprintf(errfile,
					"qhull internal warning (main): did not free %d bytes "
					"of long memory (%d pieces)\n", 
					totlong, curlong);

		return not exitcode;
	};

	template<typename T>
	polymake::Matrix<Rational> Polytope<T>::polymakePoints() const {
		polymake::Matrix<Rational> matrix(vertices.size(),
										  ambientDim + 1);
		for (int i = 0; i < matrix.rows(); ++i) {
			matrix(i, 0) = 1;
			for (int j = 1; j < matrix.cols(); ++j) {
				matrix(i, j) = vertices[i][j - 1];
			}
		}
		return matrix;
	}
	
	template<typename T>
	template<typename Solver>
	bool Polytope<T>::redundantPolymake(std::vector<bool>& redundant,
										Solver& s) const {
		polymake::Matrix<Rational> points = polymakePoints();

		polymake::Bitset keep = s.find_vertices_among_points(points).first;

		for (int i = 0; i < points.rows(); ++i) {
			redundant[i] = not keep.contains(i);
		}

		return true;
	}
	
	template<typename T>
	bool Polytope<T>::redundantCDD(std::vector<bool>& redundant) const {
		return redundantPolymake(redundant, CDDSolver);
	}

	template<typename T>	
	bool Polytope<T>::redundantLRS(std::vector<bool>& redundant) const {
		return redundantPolymake(redundant, LRSSolver);
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& strm, const Polytope<T>& p) {
		typedef typename Polytope<T>::VectorList::const_iterator VertexIterator;
		
		p.ensureReduced();
		VertexIterator v = p.vertices.begin();
		while (v != p.vertices.end()) {
			strm << *v << '\n';
			++v;
		}

		return strm;
	}
}

#endif // __POLYTOPE_POLYTOPE_HH__
