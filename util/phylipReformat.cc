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

#include "bio/formats/phylip/InputStream.hh"
#include "bio/formats/phylip/OutputStream.hh"
#include "util/options.hh"

using namespace bio::formats::phylip;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool inputSequential = false;
	bool outputSequential = false;
	size_t spaceFreq = Constants::DEFAULT_SPACE_FREQ;
	size_t lineWidth = Constants::DEFAULT_LINE_WIDTH;
	size_t blockSpacing = Constants::DEFAULT_BLOCK_SPACING;
	
	// Parse command line
	util::options::Parser parser("< phylipInput > phylipOutput",
								 "Reformat PHYLIP input");
	parser.addStoreTrueOpt(0, "in-sequential", 
						   "Input is in sequential format",
						   inputSequential);
	parser.addStoreTrueOpt(0, "out-sequential", 
						   "Output should be written in sequential format",
						   outputSequential);
	parser.addStoreOpt(0, "space-freq", 
					   "Sequence to be separated by a space every NUM "
					   "characters (no separation if NUM = 0)",
					   spaceFreq, "NUM");
	parser.addStoreOpt(0, "line-width", 
					   "Maximum length of a line to be written",
					   lineWidth, "NUM");
	parser.addStoreOpt(0, "block-spacing", 
					   "Number of blank lines between alignment blocks",
					   blockSpacing, "NUM");
	parser.parse(argv, argv + argc);

	try {
		InputStream inputStream(std::cin);
		inputStream.setSequential(inputSequential);
		OutputStream outputStream(std::cout);
		outputStream.setSequential(outputSequential);
		outputStream.setSpaceFreq(spaceFreq);
		outputStream.setLineWidth(lineWidth);
		outputStream.setBlockSpacing(blockSpacing);
		
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
