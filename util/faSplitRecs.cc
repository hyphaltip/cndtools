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

#include <string>
#include <iostream>  // for I/O
#include <algorithm> // for transform()
#include <cctype>    // for toupper()
#include <stdexcept> // for runtime_error

#include "bio/formats/fasta.hh"
#include "util/string.hh"
#include "util/stl.hh"
#include "filesystem/Path.hh"
#include "util/options.hh"
	
void writeRecord(const bio::formats::fasta::Record& rec,
				 const filesystem::Path& writeDir) {
	// Construct name of record file from first word in description
	std::string recName = util::string::validPosix(util::string::firstWord(rec.title));
	if (recName == "") {
		throw std::runtime_error("Could not write file for record with "
								 "description: " + rec.title);
	}

	// Open record file and write record
	filesystem::Path writeFilePath = writeDir / (recName + ".fa");
	std::ofstream writeFile;
	writeFilePath.openForOutput(writeFile);
	bio::formats::fasta::OutputStream fastaOutStream(writeFile);
	fastaOutStream << rec;
	std::cerr << writeFilePath.toString() << '\n';
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	std::string outDirString = ".";

	// Parse command line
	util::options::Parser parser("< fastaInput",
								 "Split FASTA file with multiple records into "
								 "multiple single-record files");
	parser.addStoreOpt('d', "out-dir",
					   "Directory in to which to output the FASTA files",
					   outDirString, "DIR");
	parser.parse(argv, argv + argc);

	try {
		filesystem::Path outDir(outDirString);
		if (not outDir.exists() or not outDir.isDirectory()) {
			throw std::runtime_error("Invalid output directory: "
									 + outDir.toString());
		}

		bio::formats::fasta::InputStream fastaStream(std::cin);
		bio::formats::fasta::Record rec;
		while (fastaStream >> rec) {
			writeRecord(rec, outDir);
		}
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
