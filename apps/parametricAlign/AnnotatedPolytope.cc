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

#include "AnnotatedPolytope.hh"

#include <Rational.h>
#include <ListMatrix.h>
#include <Vector.h>

using polymake::ListMatrix;
using polymake::Vector;

bool AnnotatedPolytope::hasSection(const std::string& title) const {
	return sections.find(title) != sections.end();
}

void AnnotatedPolytope::setSection(Section* section) {
	sections[section->title] = section;
}

void BasicSection::read(util::io::line::InputStream& stream) {
	std::string line;
	while (stream >> line and not line.empty()) {
		lines.push_back(line);
	}
}

void BasicSection::write(std::ostream& stream) {
	for (size_t i = 0; i < lines.size(); ++i) {
		stream << lines[i] << '\n';
	}
}

template<typename T>
void MatrixSection<T>::read(util::io::line::InputStream& stream) {
	ListMatrix< Vector<T> > listMatrix;
	std::string line;
	while (stream >> line and not line.empty()) {
		std::istringstream eltStream(line);
		std::vector<T> elts;
		T elt;
		while (eltStream >> elt) {
			elts.push_back(elt);
		}
		Vector<T> v(elts.size(), elts.begin());
		listMatrix /= v;
	}
	matrix = listMatrix;
}

template<typename T>
void MatrixSection<T>::write(std::ostream& stream) {
	for (int i = 0; i < matrix.rows(); ++i) {
		for (int j = 0; j < matrix.cols(); ++j) {
			if (j > 0) { stream << ' '; }
			stream << matrix(i, j);
		}
		stream << '\n';
	}
}

void SetListSection::read(util::io::line::InputStream& stream) {
	std::string line;
	while (stream >> line and not line.empty()) {
		std::istringstream intStream(line.substr(1, line.size() - 2));
		std::vector<int> indices;
		int index;
		while (intStream >> index) {
			indices.push_back(index);
		}
		std::sort(indices.begin(), indices.end());
		polymake::Set<int> s(indices.begin(), indices.end());
		sets.push_back(s);
	}
}

void SetListSection::write(std::ostream& stream) {
	for (size_t i = 0; i < sets.size(); ++i) {
		const polymake::Set<int>& s = sets[i];
		stream << '{';
		for (polymake::Set<int>::const_iterator it = s.begin(); it != s.end(); ++it) {
			if (it != s.begin()) { stream << ' '; }
			stream << *it;
		}
		stream << "}\n";
	}
}

void PeterSetListSection::read(util::io::line::InputStream& stream) {
	std::string line;
	while (stream >> line and not line.empty()) {
		std::istringstream intStream(line);
		std::vector<int> indices;
		int index;
		while (intStream >> index) {
			indices.push_back(index);
		}
		std::sort(indices.begin(), indices.end());
		polymake::Set<int> s(indices.begin(), indices.end());
		sets.push_back(s);
	}
}

void PeterSetListSection::write(std::ostream& stream) {
	for (size_t i = 0; i < sets.size(); ++i) {
		const polymake::Set<int>& s = sets[i];
		for (polymake::Set<int>::const_iterator it = s.begin(); it != s.end(); ++it) {
			if (it != s.begin()) { stream << ' '; }
			stream << *it;
		}
		stream << "\n";
	}
}

std::istream& operator>>(std::istream& stream, AnnotatedPolytope& f) {
	std::string line;
	std::string title;
	util::io::line::InputStream lineStream(stream);
	while (true) {
		while (lineStream >> line and
			   (line.empty() or line.at(0) == '_'));
		if (not lineStream) { break; }

		title = line;
		Section* section = NULL;
		if (title == "POINTS" or
			title == "VERTICES" or
			title == "FACETS" or
			title == "INEQUALITIES") {
			section = new MatrixSection<Rational>();
		} else if (title == "FACET_NORMALS" or
				   title == "REASONABLE_VERTICES") {
			section = new MatrixSection<Integer>();
		} else if (title == "FACETS_THRU_VERTICES") {
			section = new SetListSection();
		} else if (title == "POINTS_IN_FACETS") {
			section = new PeterSetListSection();
		} else {
			section = new BasicSection();
		}
		section->title = title;
		section->read(lineStream);
		f.setSection(section);
	}
	return stream;
}

std::ostream& operator<<(std::ostream& stream, AnnotatedPolytope& f) {
	for (AnnotatedPolytope::SectionMap::const_iterator it = f.sections.begin();
		 it != f.sections.end(); ++it) {
		Section* section = it->second;
		stream << section->title << '\n';
		section->write(stream);
		stream << '\n';
	}
	return stream;
}
