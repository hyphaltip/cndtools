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
#include "bio/formats/fasta.hh"
#include "util/options.hh"

using namespace bio::formats;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool outputSequential = false;
	size_t spaceFreq = phylip::Constants::DEFAULT_SPACE_FREQ;
	size_t lineWidth = phylip::Constants::DEFAULT_LINE_WIDTH;
	size_t blockSpacing = phylip::Constants::DEFAULT_BLOCK_SPACING;
	
	// Parse command line
	util::options::Parser parser("< fastaInput > phylipOutput",
								 "Convert multi-FASTA alignments to PHYLIP alignments");
	parser.addStoreTrueOpt('s', "out-sequential", 
						   "PHYLIP output should be written in sequential format",
						   outputSequential);
	parser.addStoreOpt('f', "space-freq", 
					   "Sequence to be separated by a space every NUM "
					   "characters in PHYLIP output (no separation if NUM = 0)",
					   spaceFreq, "NUM");
	parser.addStoreOpt('w', "line-width", 
					   "Maximum length of a line to be written in PHYLIP output",
					   lineWidth, "NUM");
	parser.addStoreOpt('b', "block-spacing", 
					   "Number of blank lines between PHYLIP alignment blocks",
					   blockSpacing, "NUM");
	parser.parse(argv, argv + argc);

	try {
		fasta::InputStream inputStream(std::cin);
		phylip::OutputStream outputStream(std::cout);
		outputStream.setSequential(outputSequential);
		outputStream.setSpaceFreq(spaceFreq);
		outputStream.setLineWidth(lineWidth);
		outputStream.setBlockSpacing(blockSpacing);
		
		// Step through records, possibly cleaning, and then output
		fasta::Record inputRec;
		phylip::Record outputRec;
		while (inputStream >> inputRec) {
			outputRec.names.push_back(inputRec.title);
			outputRec.sequences.push_back(inputRec.sequence);
		}
		outputStream << outputRec;
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
