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

#ifndef __ANNOTATED_POLYTOPE_HH__
#define __ANNOTATED_POLYTOPE_HH__

#include "util/io/line/InputStream.hh"
#include "util/stl.hh"

#include <Matrix.h>
#include <Set.h>

class Section;

class AnnotatedPolytope {
public:
	template<typename T>
	void getSection(const std::string& title, T*& section) const;

	bool hasSection(const std::string& title) const;
	void setSection(Section* section);

	friend std::ostream& operator<<(std::ostream& stream, AnnotatedPolytope& f);
		
private:
	typedef util::stl::hash_map<std::string, Section*> SectionMap;
	SectionMap sections;
};

std::istream& operator>>(std::istream& stream, AnnotatedPolytope& f);
std::ostream& operator<<(std::ostream& stream, AnnotatedPolytope& f);

class Section {
public:
	virtual ~Section() {}
	Section(const std::string& title = "") : title(title) {}
	virtual void read(util::io::line::InputStream& stream) = 0;
	virtual void write(std::ostream& stream) = 0;
	std::string title;
};

class BasicSection : public Section {
public:
	std::vector<std::string> lines;
	void read(util::io::line::InputStream& stream);
	void write(std::ostream& stream);
};

template<typename T>
class MatrixSection : public Section {
public:
	polymake::Matrix<T> matrix;
	void read(util::io::line::InputStream& stream);
	void write(std::ostream& stream);
};

class SetListSection : public Section {
public:
	std::vector< polymake::Set<int> > sets;
	virtual ~SetListSection() {}
	virtual void read(util::io::line::InputStream& stream);
	virtual void write(std::ostream& stream);
};

class PeterSetListSection : public SetListSection {
public:
	void read(util::io::line::InputStream& stream);
	void write(std::ostream& stream);
};

template<typename T>
void AnnotatedPolytope::getSection(const std::string& title,
								   T*& section) const {
	if (not hasSection(title)) {
		throw std::runtime_error("Polytope does not have section " + title);
	}
	section = dynamic_cast<T*>(sections.find(title)->second);
	if (section == NULL) {
		throw std::runtime_error("Section " + title + " not of specified type");
	}
}

#endif // __ANNOTATED_POLYTOPE_HH__
