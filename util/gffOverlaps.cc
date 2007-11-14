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
#include <vector>
#include <set>
#include <algorithm>

#include "boost/function_output_iterator.hpp"

#include "bio/gff/GFFRecord.hh"
#include "bio/gff/GFFInputStream.hh"
#include "bio/genome/Coord.hh"
#include "util/interval.hh"
#include "util/options.hh"
#include "filesystem.hh"

using namespace filesystem;
using bio::gff::GFFRecord;
using bio::gff::GFFInputStream;

template<typename T>
struct DereferencedLessThan {
	bool operator()(const T* x, const T* y) const { return *x < *y; }
};

namespace util { namespace interval {
	template<>
	struct interval_traits<GFFRecord*> {
		typedef bio::genome::Coord coord_type;
		typedef bio::genome::Distance difference_type;
		static coord_type start(const GFFRecord* i) {
			return bio::genome::Coord(i->getSeqname(), i->getStart() - 1);
		}
		static coord_type end(const GFFRecord* i) {
			return bio::genome::Coord(i->getSeqname(), i->getEnd());
		}
	};
} }

void readGFFRecords(const std::string& filename,
					std::vector<GFFRecord*>& recs) {
	InputFileStream file(filename);
	GFFInputStream stream(file);
	GFFRecord r;
	while (stream >> r) {
		recs.push_back(new GFFRecord(r));
	}
}

void writeGFFRecords(const std::string& filename,
					 const std::vector<GFFRecord*>& recs) {
	OutputFileStream file(filename);
	typedef std::vector<GFFRecord*>::const_iterator Iterator;
	for (Iterator it = recs.begin(); it != recs.end(); ++it) {
		file << **it;
	}
}

void sortGFFRecords(std::vector<GFFRecord*>& recs) {
	std::sort(recs.begin(), recs.end(), DereferencedLessThan<GFFRecord>());
}

// Functor for receiving the results of the overlap analysis of the
// GFF intervals.
struct OverlapRecorder {
	OverlapRecorder(std::set<GFFRecord*>& overlapped1,
					std::set<GFFRecord*>& overlapped2)
		: overlapped1(overlapped1), overlapped2(overlapped2) {
	}

 	void operator()(const std::pair<GFFRecord*, GFFRecord*>& overlap) {
		overlapped1.insert(overlap.first);
		overlapped2.insert(overlap.second);
	}
	
	std::set<GFFRecord*>& overlapped1;
	std::set<GFFRecord*>& overlapped2;
};

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Options and defaults
	std::string gffInFilename1;
	std::string gffInFilename2;
	std::string gffOutFilename1;
	std::string gffOutFilename2;
	
	util::options::Parser parser("", "");
	parser.addStoreArg("gffFile1", "First GFF File", gffInFilename1);
	parser.addStoreArg("gffFile2", "Second GFF File", gffInFilename2);
	parser.addStoreArg("outGFFFile1", "First output GFF File", gffOutFilename1);
	parser.addStoreArg("outGFFFile2", "Second output GFF File", gffOutFilename2);
	parser.parse(argv, argv + argc);

	try {
		std::vector<GFFRecord*> gffRecs1, gffRecs2;
		readGFFRecords(gffInFilename1, gffRecs1);
		readGFFRecords(gffInFilename2, gffRecs2);

		sortGFFRecords(gffRecs1);
		sortGFFRecords(gffRecs2);

		std::set<GFFRecord*> overlapped1, overlapped2;
		util::interval::overlaps(gffRecs1.begin(), gffRecs1.end(),
								 gffRecs2.begin(), gffRecs2.end(),
								 boost::make_function_output_iterator(OverlapRecorder(overlapped1, overlapped2)));


		std::vector<GFFRecord*> overlappedList1(overlapped1.begin(),
												overlapped1.end());
		std::vector<GFFRecord*> overlappedList2(overlapped2.begin(),
												overlapped2.end());

		sortGFFRecords(overlappedList1);
		sortGFFRecords(overlappedList2);
		
		writeGFFRecords(gffOutFilename1, overlappedList1);
		writeGFFRecords(gffOutFilename2, overlappedList2);

	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what();
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
