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

#include "bio/alignment/BasicNamedMultipleAlignment.hh"
#include "bio/formats/fasta.hh"
#include "util/options.hh"
using namespace bio::alignment;
using namespace bio::formats;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);
	
	// Initialize options to defaults
	size_t colStart;
	size_t colEnd;
	
	// Parse command line
	util::options::Parser parser("< mfaInput", "Slice multiple alignments");
	parser.addStoreArg("colStart", "first column to include in slice",
					   colStart);
	parser.addStoreArg("colEnd", "last column to include in slice",
					   colEnd);
	parser.parse(argv, argv + argc);

	try {
		// Read in multiple alignment
		fasta::InputStream fastaStream(std::cin);
		BasicNamedMultipleAlignment inAlignment;
		fastaStream >> inAlignment;
		
		// Check for valid columns
		if (colEnd >= inAlignment.getNumCols()) {
			throw std::runtime_error("end column is greater than number of "
									 "columns in alignment");
		} else if (colStart > colEnd) {
			throw std::runtime_error("start column is greater than end column");
		}

		// Construct alignment slice
		BasicNamedMultipleAlignment outAlignment;
		for (size_t i = 0; i < inAlignment.getNumSeqs(); ++i) {
			outAlignment.addSeq(inAlignment.getSubstring(i,
														 colStart,
														 colEnd),
								inAlignment.getName(i));
		}

		// Output alignment slice
		fasta::OutputStream outStream(std::cout);
		outStream << outAlignment;
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
