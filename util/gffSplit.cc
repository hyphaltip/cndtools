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
#include <fstream>

#include "bio/gff/GFFRecord.hh"
#include "bio/gff/GFFInputStream.hh"
#include "util/parser.hh"
#include "util/stl.hh"
#include "util/options.hh"
#include "filesystem/Path.hh"

int main(int argc, const char** argv) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);
	
	// Values to be (optionally) specified on the command line
	bool bySeqname = false;
	bool bySource = false;
	bool byFeature = false;
	std::string outdir = ".";
	
	// Set up options
	util::options::Parser parser("< gffInput", "");
	parser.addStoreTrueOpt('c', "by-seqname",
						   "split into separate files by seqname field",
						   bySeqname);
	parser.addStoreTrueOpt('s', "by-source",
						   "split into separate files by source field (default)",
						   bySource);
	parser.addStoreTrueOpt('f', "by-feature",
						   "split into separate files by feature field",
						   byFeature);
	parser.addStoreOpt('o', "outdir", 
					   "directory to output split gff files",
					   outdir, "DIR");
	parser.parse(argv, argv + argc);

	// Use source as default field on which to split
	if (!(bySeqname || bySource || byFeature)) {
		bySource = 1;
	}

	// Check that at only one field is specified
	if ((bySeqname + bySource + byFeature) > 1) {
		throw std::runtime_error("May only specify one field on which to split");
	}
	
	// Form output directory path and check for its existence
	filesystem::Path outDir(outdir);
	if (not outDir.exists()) {
		std::cerr << "Error: Invalid output directory: " << outdir << '\n';
		return EXIT_FAILURE;
	}
	
	// Make a map from field values to files for records with that field value
	typedef util::stl::hash_map<std::string, std::ofstream*> FileMap;
	FileMap files;
	
	try {
		// Construct GFF stream for fast reading
		bio::gff::GFFInputStream gffStream(std::cin);

		bio::gff::GFFRecord rec;
		while (gffStream >> rec) {

			// Extract the appropriate field
			std::string field;
			if (bySeqname) {
				field = rec.getSeqname();
			} else if (bySource) {
				field = rec.getSource();
			} else { // byFeature
				field = rec.getFeature();
			}

			// Output records that do not have the field on standard out
			if (field == ".") {
				std::cout << rec;
				continue;
			}

			// Get the file for this field.  Create a new file if no
			// file has been opened for this field value
			std::ofstream* outFile;			
			FileMap::iterator pos = files.find(field);
			if (pos == files.end()) {
				outFile = new std::ofstream();
				(outDir / (field + ".gff")).openForOutput(*outFile);
				pos = files.insert(std::make_pair(field, outFile)).first;
			} else {
				outFile = pos->second;
			}

			// Send this record to the appropriate file
			(*outFile) << rec;
		}
	} catch (util::parser::FormatError& e) {
		std::cerr << e.getProblem() << '\n' << e.getLine() << '\n';
	}

	// Close all of the opened files
	for (FileMap::iterator pos = files.begin(); pos != files.end(); ++pos) {
		delete pos->second;
	}
}
