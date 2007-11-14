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
#include <algorithm>

#include "util/options.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	size_t n;

	// Parse command line
	util::options::Parser parser("", "");
	parser.addStoreArg("N", "Length of the permutation", n);
	parser.parse(argv, argv + argc);

	std::vector<size_t> p(n);
	for (size_t i = 0; i < n; ++i) { p[i] = i; }
	std::random_shuffle(p.begin(), p.end());
	for (size_t i = 0; i < n; ++i) { std::cout << p[i] << '\n'; }
	
	return EXIT_SUCCESS;
}
