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

#include "util/stl.hh"
#include "util/options.hh"

const std::string USAGE = "";

const std::string DESCRIPTION =
	"Generates N^D grid points in the unit hypercube of dimension D";

void print_grid(int n, size_t d,
				size_t curr_dim, std::vector<double>& curr_point,
				std::ostream& stream) {
	for (int i = 0; i < n; ++i) {
		curr_point[curr_dim] = (static_cast<double>(i) + 0.5) / n;
		if (curr_dim == d - 1) {
			util::stl::print_elements(stream, curr_point, "\t");
		} else {
			print_grid(n, d, curr_dim + 1, curr_point, stream);
		}
	}
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	size_t n = 10;
	size_t d = 2;
	
	// Set up option parser
	util::options::Parser parser(USAGE, DESCRIPTION);
	parser.addStoreArg("n", "Number of grid values in each dimension", n);
	parser.addStoreOpt('d', "dimension",
					   "Dimension of hypercube in which to generate points",
					   d, "INTEGER");
	parser.parse(argv, argv + argc);

	try {

		std::vector<double> point(d);
		print_grid(n, d, 0, point, std::cout);
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
