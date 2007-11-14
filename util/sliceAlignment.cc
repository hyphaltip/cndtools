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

#include "bio/alignment/AlignmentSlicer.hh"
#include "bio/genome/BasicInterval.hh"
#include "bio/formats/fasta.hh"
#include "util/options.hh"
#include "util/string.hh"
#include "util/io/line/InputStream.hh"
#include "filesystem.hh"
using namespace bio::formats;
using util::string::toString;
using bio::alignment::AlignmentSlicer;
using namespace util::io;

const std::string description =
"Extracts alignments out of a whole genome multiple alignment.  The multiple "
"alignment must be in the form of an orthology map combined with alignments "
"for each segment in the orthology map.  The specified alignment "
"directory must contain an orthology map file named \"map\" as well as "
"a file \"genomes\" "
"that indicates the order of the genomes in the map file. The alignment "
"directory must contain a numbered subdirectory for each segment in the "
"orthology map.  The segment subdirectories must each contain a file named "
"\"mavid.mfa\" that has a multiple alignment of that segment in multi-FASTA "
"format.\n\n"
"A single interval to extract from the alignment may be specified as command "
"line arguments.  If multiple intervals are to be extracted, they may be "
"given on the standard input, with one interval specified per line. "
"Intervals must be of the form CHROM START END STRAND and should be "
"zero-based half-open [START, END). "
"The output consists of multiple alignments for each input interval in "
"multi-FASTA format.  Every multiple alignment will have the same number of "
"records (the number of genomes in the map), allowing the user to separate "
"the different alignments.";

typedef std::vector<std::string> Coords;

util::string::Converter<bio::genome::Position> to_position;

std::string to_string(const Coords& coords) {
	return util::string::join(coords.begin(), coords.end(), " ");
}

bio::genome::BasicInterval to_interval(const Coords& coords) {
	if (coords.size() != 4 or coords[3].size() != 1) {
		throw std::runtime_error("Invalid interval: " + to_string(coords));
	}
	bio::genome::BasicInterval interval(coords[0],
										to_position(coords[1]),
										to_position(coords[2]),
										coords[3][0]);
	if (interval.getStart() >= interval.getEnd()) {
		throw std::runtime_error("Invalid interval: " + to_string(coords));
	}
	return interval;
}

line::InputStream& operator>>(line::InputStream& stream, Coords& coords) {
	static std::string line;
	if (stream >> line) {
		coords.clear();
		util::string::split(line, std::back_inserter(coords));
	}
	return stream;
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	char noMapChar = '-';
	std::string alignDirname;
	std::string sourceGenome;
	Coords coords;

	util::options::Parser parser("[chrom start end strand] [< intervals]",
								 description);
	parser.addStoreOpt(0, "no-map-char",
					   "character to use in the alignment for intervals for "
					   "which there is no orthologous segment identified in a "
					   "given genome",
					   noMapChar, "CHAR");
	parser.addStoreArg("align_dir", "directory containing alignments",
					   alignDirname);
	parser.addStoreArg("genome",
					   "Name of reference genome for which coordinates "
					   "are specified", sourceGenome);
	parser.addAppendArg("", "", coords, 0, 4);
	parser.parse(argv, argv + argc);

	try {
		AlignmentSlicer slicer(alignDirname, sourceGenome, noMapChar);
		
		fasta::OutputStream fastaOutStream(std::cout);

		if (coords.empty()) {
			util::io::line::InputStream lineStream(std::cin);
			while (lineStream >> coords) {
				fastaOutStream << slicer.getSlice(to_interval(coords));
			}
		} else {
			fastaOutStream << slicer.getSlice(to_interval(coords));
		}
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
