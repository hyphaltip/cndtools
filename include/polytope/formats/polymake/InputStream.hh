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

#ifndef __POLYTOPE_FORMATS_POLYMAKE_INPUTSTREAM_HH__
#define __POLYTOPE_FORMATS_POLYMAKE_INPUTSTREAM_HH__

#include <iosfwd>
#include <sstream>

#include "util/io/line/InputStream.hh"
#include "polytope/Polytope.hh"

namespace polytope { namespace formats { namespace polymake {

	class InputStream {
	public:
		InputStream(std::istream& strm,
					const size_t bufferSize=IO_BUFFER_SIZE);

		template<typename T>
		InputStream& operator>>(Polytope<T>& p);

		operator bool();
		bool operator!();
		
	private:
		void readToNextSection();
		
		void skipSection();
		
		template<typename T>
		Vector<T> readVector(std::istream& strm);
		
		template<typename T>
		void readVertices(typename Polytope<T>::VectorList& vl);
		
		template<typename T>
		void readVertexCounts(typename Polytope<T>::CountList& cl);
		
		util::io::line::InputStream strm;
		std::string line;
	};

	inline InputStream::InputStream(std::istream& strm,
									const size_t bufferSize) :
		strm(strm, bufferSize) {
	}
	
	inline InputStream::operator bool() {
		return strm;
	}
	
	inline bool InputStream::operator!() {
		return !strm;
	}

	void InputStream::readToNextSection() {
		while (strm and strm >> line and line.empty());
	}

	void InputStream::skipSection() {
		while (strm and strm >> line and not line.empty());
	}
	
	template<typename T>
	Vector<T> InputStream::readVector(std::istream& strm) {
		std::vector<T> coords;
		T coord;
		while(strm >> coord) { coords.push_back(coord); }
		Vector<T> v(coords.size() - 1);
		for (size_t i = 0; i < v.size(); ++i) {
			v[i] = coords[i + 1];
		}
		return v;
	}
	
	template<typename T>
	void InputStream::readVertices(typename Polytope<T>::VectorList& vl) {
		while (strm and strm >> line and not line.empty()) {
			std::istringstream lineStream(line);
			vl.push_back(readVector<T>(lineStream));
		}
	}

// 	template<typename T>
// 	void InputStream::readVertexCounts(typename Polytope<T>::CountList& cl) {
// 		while (strm and strm >> line and not line.empty()) {
// 			std::istringstream lineStream(line);
// 			typename Polytope<T>::Count c;
// 			lineStream >> c;
// 			cl.push_back(c);
// 		}
// 	}
	
	template<typename T>
	InputStream& InputStream::operator>>(Polytope<T>& p) {
		typename Polytope<T>::VectorList vl;
// 		typename Polytope<T>::CountList cl;

		while (strm) {
			readToNextSection();
			if (line == "VERTICES") { readVertices<T>(vl); }
// 			else if (line == "VERTEX_COUNTS") { readVertexCounts<T>(cl); }
			else { skipSection(); }
		}

// 		if (cl.empty()) {
			p = Polytope<T>(vl, true);
// 		} else {
// 			p = Polytope<T>(vl, cl, true);
// 		}
		
		return *this;
	}

	
} } }

#endif // __POLYTOPE_FORMATS_POLYMAKE_INPUTSTREAM_HH__
