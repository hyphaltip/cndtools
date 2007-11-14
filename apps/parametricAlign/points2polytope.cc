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
#include <sstream>

#include "polytope/Polytope.hh"
#include "polytope/formats/polymake/OutputStream.hh"
#include "util/options.hh"
#include "util/io/line/InputStream.hh"
using namespace polytope;
using namespace util::io;

typedef int T;
typedef Polytope<T> Poly;
typedef Vector<T> Point;

void readPoints(std::istream& stream, Poly::VectorList& points) {
	line::InputStream lineInputStream(stream);
	std::string line;
	std::vector<T> coords;
	while (lineInputStream >> line) {
		coords.clear();
		std::istringstream coordStream(line);
		T coord;
		while (coordStream >> coord) {
			coords.push_back(coord);
		}
		Point point(coords.size());
		for (size_t i = 0; i < coords.size(); ++i) {
			point[i] = coords[i];
		}
		points.push_back(point);
	}
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments

	// Set up option parser
	util::options::Parser parser("", "Converts a list of points to a polytope");
	parser.parse(argv, argv + argc);

	try {
		// Read in points
		Poly::VectorList points;
		readPoints(std::cin, points);

		// Make polytope from points
		Poly polytope(points);

		// Output polytope
		formats::polymake::OutputStream polymakeOutputStream(std::cout);
		polymakeOutputStream << polytope;

	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
