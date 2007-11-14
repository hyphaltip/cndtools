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

#include "bio/formats/fasta.hh"
#include "util/stl.hh"
#include "util/options.hh"
#include "filesystem.hh"
using util::stl::hash_set;
using namespace bio::formats;
using namespace filesystem;

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	std::vector<std::string> titlesList;
	std::string titleFilename;
	
	util::options::Parser parser("< fastaInput",
								 "Outputs records in the input with titles "
								 "given as arguments");
	parser.addStoreOpt('f', "file", "file containing titles", titleFilename);
	parser.addAppendArg("title", "title of record to keep in output",
						titlesList);
	parser.parse(argv, argv + argc);

	try {
		if (not titleFilename.empty()) {
			InputFileStream titleFile(titleFilename);
			std::string title;
			while (std::getline(titleFile, title)) {
				titlesList.push_back(title);
			}
		}
	
		hash_set<std::string> titles(titlesList.begin(), titlesList.end());
	
		// Construct FASTA stream for fast reading
		fasta::InputStream fastaInputStream(std::cin);
		fasta::OutputStream fastaOutputStream(std::cout);

		// Step through records, output appropriate records
		size_t numFound = 0;
		fasta::Record rec;
		while (numFound != titles.size() and fastaInputStream >> rec) {
			if (titles.find(rec.title) != titles.end()) {
				fastaOutputStream << rec;
				++numFound;
			}
		}
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
