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
#include <algorithm>
#include <map>
#include <vector>

#include "bio/formats/fasta.hh"
#include "bio/formats/repeatmasker/InputStream.hh"
#include "util/string.hh"
#include "util/options.hh"
#include "filesystem.hh"

using namespace filesystem;
using namespace bio;
using namespace bio::formats;

struct LowerCaser {
	char operator()(const char c) const { return std::tolower(c); }
};

class RepeatMasker {
public:
	RepeatMasker(repeatmasker::InputStream& strm);
	void mask(const std::string& name, std::string& seq) const;

private:
	void mask(std::string& seq, const genome::Interval& i) const;
	
	typedef std::vector<repeatmasker::Record> RecordList;
	typedef std::map<std::string, RecordList> RecordListMap;
	RecordListMap intervals;
};

RepeatMasker::RepeatMasker(repeatmasker::InputStream& strm) {
	repeatmasker::Record rec;
	while (strm >> rec) {
		intervals[rec.getChrom()].push_back(rec);
	}
}

void RepeatMasker::mask(std::string& seq, const genome::Interval& i) const {
	std::transform(seq.begin() + i.getStart(),
				   seq.begin() + i.getEnd(),
				   seq.begin() + i.getStart(),
				   LowerCaser());
}

void RepeatMasker::mask(const std::string& name, std::string& seq) const {
	RecordListMap::const_iterator entry = intervals.find(name);
	if (entry == intervals.end()) { return; }
	const RecordList& intervalList = entry->second;
	for (RecordList::const_iterator it = intervalList.begin();
		 it != intervalList.end(); ++it) {
		mask(seq, *it);
	}
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults
	std::string rmaskFilename;

	// Parse command line
	util::options::Parser parser("< fastaInput",
								 "Softmask sequence according to RepeatMasker"
								 "tables");
	parser.addStoreArg("rmaskFile", "RepeatMasker table file", rmaskFilename);
	parser.parse(argv, argv + argc);
	
	try {
		// Open RepeatMasker table file
		InputFileStream rmaskFile(rmaskFilename);
		repeatmasker::InputStream rmaskInStream(rmaskFile);
		RepeatMasker masker(rmaskInStream);

		// Construct FASTA stream for fast reading
		fasta::InputStream fastaInStream(std::cin);
		fasta::OutputStream fastaOutStream(std::cout);
		
		// Step through records, possibly cleaning, and then output
		fasta::Record rec;
		while (fastaInStream >> rec) {
			masker.mask(rec.title, rec.sequence);
			fastaOutStream << rec;
		}
		
	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
		
	return EXIT_SUCCESS;
}
