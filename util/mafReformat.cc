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

#include "bio/formats/maf/InputStream.hh"
#include "bio/formats/maf/OutputStream.hh"
#include "util/options.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults

	// Parse command line
	util::options::Parser parser("< mafInput", "Reformat MAF input");
	parser.parse(argv, argv + argc);

	bio::formats::maf::InputStream inputStream(std::cin);
	bio::formats::maf::OutputStream outputStream(std::cout,
												 inputStream.getHeader());

	// Step through records, possibly cleaning, and then output
	bio::formats::maf::Record rec;
    while (inputStream >> rec) {
		outputStream << rec;
	}

	return EXIT_SUCCESS;
}
