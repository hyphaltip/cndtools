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

#include "bio/formats/fasta/InputStream.hh"
#include "bio/sdb.hh"
#include "util/options.hh"
#include "filesystem/Path.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options and arguments
	bool append = false;
	bool compressed = false;
	std::string sdbFilename;

	// Parse command line
	util::options::Parser parser("< fastaInput", "");
	parser.addStoreTrueOpt('a', "append",
						   "append records to an existing database",
						   append);
	parser.addStoreTrueOpt('c', "compressed",
						   "use nib compression (DNA only)",
						   compressed);
	parser.addStoreArg("sdbFile",
					   "SDB file to create or append to",
					   sdbFilename);
	parser.parse(argv, argv + argc);

	try {
		bio::SDB::DB db;
		db.open(sdbFilename, false, not append);
	
		// Construct FASTA stream for fast reading
		bio::formats::fasta::InputStream fastaStream(std::cin);
		
		// Read FASTA records into database
		bio::formats::fasta::Record rec;
		while (fastaStream >> rec) {
			db.putRec(rec.title, rec.sequence, compressed);
		}
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
