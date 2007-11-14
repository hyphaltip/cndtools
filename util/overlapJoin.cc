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

#include "util/string.hh"
#include "util/options.hh"
#include "bio/genome/Interval.hh"
#include "util/io/line/InputStream.hh"
#include "util/interval.hh"
#include "filesystem.hh"
#include "boost/function_output_iterator.hpp"
#include "boost/shared_ptr.hpp"

using namespace bio;
using namespace filesystem;

struct IntervalLine {
	std::string line;
	genome::Coord start;
	genome::Coord end;
};

struct IntervalLineTraits {
	typedef genome::Coord coord_type;
	static coord_type start(const IntervalLine& x) { return x.start; }
	static coord_type end(const IntervalLine& x) { return x.end; }
};

class IntervalLineIterator {
public:
	typedef IntervalLine value_type;
		
	IntervalLineIterator() :
		lineStream(),
		valid(false) {
	}

	IntervalLineIterator(std::istream& stream,
						 const std::vector<size_t>& indices) :
		lineStream(new util::io::line::InputStream(stream)),
		contigIndex(indices[0]),
		startIndex(indices[1]),
		endIndex(indices[2]),
		valid(true) {
		++(*this);
	}

	IntervalLine& operator*() {
		return intervalLine;
	}

	IntervalLineIterator& operator++() {
		assert(lineStream != NULL);
		(*lineStream) >> line;
		valid = (*lineStream);
		if (valid) {
			intervalLine.line = line;
			fields.clear();
			util::string::split(line, std::back_inserter(fields), "\t");
			util::string::Converter<genome::Position> toGenomePosition;
			intervalLine.start = genome::Coord(fields[contigIndex],
											   toGenomePosition(fields[startIndex]));
			intervalLine.end = genome::Coord(fields[contigIndex],
											 toGenomePosition(fields[endIndex]));
		}
		return *this;
	}

	bool operator==(const IntervalLineIterator& x) const {
		return ((not valid and not x.valid) or
				(valid and x.valid and lineStream == x.lineStream));
	}

	bool operator!=(const IntervalLineIterator& x) const {
		return not (*this == x);
	}
	
private:
	typedef boost::shared_ptr<util::io::line::InputStream> LineStreamPtr;
	
	LineStreamPtr lineStream;
	size_t contigIndex;
	size_t startIndex;
	size_t endIndex;
	IntervalLine intervalLine;
	bool valid;
	std::string line;
	std::vector<std::string> fields;
};

std::vector<size_t> parse_field_indices(const std::string& s) {
	std::vector<std::string> string_indices;
	util::string::split(s, std::back_inserter(string_indices), ",");
	std::vector<size_t> indices;
	util::string::Converter<size_t> toNum;
	for (size_t i = 0; i < string_indices.size(); ++i) {
		indices.push_back(toNum(string_indices[i]) - 1);
	}
	return indices;
}

struct OverlapRecorder {
	void operator()(const std::pair<IntervalLine, IntervalLine>& p) const {
		std::cout << p.first.line << '\t' << p.second.line << '\n';
	}
};

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::string fields1 = "1,2,3";
	std::string fields2 = "1,2,3";
	std::string filename1;
	std::string filename2;
	
	// Parse command line
	util::options::Parser parser("", "");
	parser.addStoreOpt('1', "",
					   "Interval fields for the first file",
					   fields1);
	parser.addStoreOpt('2', "",
					   "Interval fields for the second file",
					   fields2);
	parser.addStoreArg("file1", "File 1", filename1);
	parser.addStoreArg("file1", "File 2", filename2);
	parser.parse(argv, argv + argc);

	try {
		std::vector<size_t> indices1(parse_field_indices(fields1));
		std::vector<size_t> indices2(parse_field_indices(fields2));

		InputFileStream file1(filename1);
		InputFileStream file2(filename2);
		
		IntervalLineIterator lineIterator1(file1, indices1);
		IntervalLineIterator lineIterator2(file2, indices2);
		IntervalLineIterator endLineIterator;
		
		util::interval::overlaps(lineIterator1, endLineIterator,
								 lineIterator2, endLineIterator,
								 boost::make_function_output_iterator(OverlapRecorder()),
								 IntervalLineTraits(),
								 IntervalLineTraits());
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
