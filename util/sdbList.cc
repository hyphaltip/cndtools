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
#include <string>
#include <stdexcept>

#include "bio/sdb.hh"
#include "util/options.hh"

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	bool number_title = false;
	bool length = false;
	std::string dbFilename;

	util::options::Parser parser("", "");
	parser.addStoreTrueOpt('n', "number",
						   "output record numbers in addition to record titles",
						   number_title);
	parser.addStoreTrueOpt('l', "length",
						   "output sequence lengths",
						   length);
	parser.addStoreArg("dbFile", "", dbFilename);
	parser.parse(argv, argv + argc);

	bio::SDB::DB db;
	
	try {
		db.open(dbFilename);
	} catch (std::runtime_error& err) {
		std::cerr << "Error: Could not open database: " << err.what() << '\n';
		exit(EXIT_FAILURE);
	}

	for (bio::SDB::DB::Iterator it = db.begin(); it != db.end(); ++it) {
		std::cout << it->getTitle();
		if (number_title) {
			std::cout << '\t' << it->getRecNum();
		}
		if (length) {
			std::cout << '\t' << it->getLength();
		}
		std::cout << '\n';
	}

	return EXIT_SUCCESS;
}
