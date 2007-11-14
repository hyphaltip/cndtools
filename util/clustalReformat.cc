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

#include "bio/formats/clustal.hh"
#include "util/options.hh"

using namespace bio::formats::clustal;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	size_t lineWidth = Constants::DEFAULT_LINE_WIDTH;
	size_t blockSpacing = Constants::DEFAULT_BLOCK_SPACING;
	size_t minNameSeqSpacing = Constants::DEFAULT_MIN_NAME_SEQ_SPACING;
	size_t minGutterLen = Constants::DEFAULT_MIN_GUTTER_LEN;
	bool seqnos = false;

	// Parse command line
	util::options::Parser parser("< clustalInput > clustalOutput",
								 "Reformat CLUSTAL input");
	parser.addStoreTrueOpt(0, "seqnos", 
						   "Output sequence positions at end of each line",
						   seqnos);
	parser.addStoreOpt(0, "line-width", 
					   "Maximum length of a line to be written",
					   lineWidth, "NUM");
	parser.addStoreOpt(0, "block-spacing", 
					   "Number of blank lines between alignment blocks",
					   blockSpacing, "NUM");
	parser.addStoreOpt(0, "min-name-seq-spacing",
					   "Minimum spacing between sequence name and sequence",
					   minNameSeqSpacing, "NUM");
	parser.addStoreOpt(0, "min-gutter-len", 
					   "Minimum width of the name gutter",
					   minGutterLen, "NUM");
	parser.parse(argv, argv + argc);

	try {
		// Create streams
		InputStream inputStream(std::cin);
		OutputStream outputStream(std::cout);

		// Set output options
		outputStream.setLineWidth(lineWidth);
		outputStream.setBlockSpacing(blockSpacing);
		outputStream.setMinNameSeqSpacing(minNameSeqSpacing);
		outputStream.setMinGutterLen(minGutterLen);
		outputStream.setSeqNos(seqnos);
		
		// Step through records, possibly cleaning, and then output
		Record rec;
		while (inputStream >> rec) {
			outputStream << rec;
		}
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
