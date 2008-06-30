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
#include <string>
#include <limits>
#include <ctime>

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int.hpp"
#include "boost/random/variate_generator.hpp"

#include "util/options.hh"

typedef boost::mt19937 base_generator_type;
typedef boost::uniform_int<int> distribution_type;
typedef boost::variate_generator<base_generator_type&,
								 distribution_type> variate_generator_type;

const std::string USAGE = "";

const std::string DESCRIPTION =
	"Generates N random integers in the range [0, 2^31 - 1]";

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Set up options and arguments
	size_t n;
	base_generator_type::result_type seed = std::time(0);
	
	// Set up option parser
	util::options::Parser parser(USAGE, DESCRIPTION);
	parser.addStoreArg("n", "Number of random integers to generate", n);
	parser.addStoreOpt('s', "seed", "Seed for random number generator",
					   seed, "INTEGER");
	parser.parse(argv, argv + argc);

	try {

		base_generator_type generator(seed);
		distribution_type uni_dist(0, std::numeric_limits<int>::max());
		variate_generator_type uni(generator, uni_dist);

		for (size_t i = 0; i < n; ++i) {
			std::cout << uni() << '\n';
		}
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
