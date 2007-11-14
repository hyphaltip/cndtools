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

#include "util/options.hh"
#include "bio/genome/BasicInterval.hh"
using namespace bio::genome;

typedef std::vector<BasicInterval> IntervalList;

void readIntervals(std::istream& stream,
				   IntervalList& fIntervals,
				   IntervalList& rIntervals) {
	std::string chrom;
	Position start, end;
	Strand strand;
	while (stream >> chrom >> start >> end >> strand) {
		BasicInterval i(chrom, start, end, strand);
		if (strand.isForward()) {
			fIntervals.push_back(i);
		} else {
			rIntervals.push_back(i);
		}
	}
}

void writeIntervals(std::ostream& stream,
					IntervalList& intervals) {
	for (size_t i = 0; i < intervals.size(); ++i) {
		BasicInterval& in = intervals[i];
		stream << in.getChrom() << '\t'
			   << in.getStart() << '\t'
			   << in.getEnd() << '\t'
			   << in.getStrand() << '\n';
	}
}

void unionIntervals(IntervalList& intervals) {
	sort(intervals.begin(), intervals.end());
	
	IntervalList::iterator prev(intervals.begin()), next(intervals.begin());
	while (next != intervals.end()) {
		if (prev->overlaps(*next)) {
			*prev |= *next;
		} else {
			++prev;
			*prev = *next;
		}
		++next;
	}

	if (prev != intervals.end()) {
		++prev;
		intervals.erase(prev, intervals.end());
	}
}

void mergeIntervals(IntervalList& fIntervals,
					IntervalList& rIntervals,
					IntervalList& intervals) {
	std::merge(fIntervals.begin(), fIntervals.end(),
			   rIntervals.begin(), rIntervals.end(),
			   std::back_inserter(intervals));
}

int main(int argc, const char* argv[]) {
	// Increase speed of input/output to standard streams
	std::ios::sync_with_stdio(false);

	// Initialize options to defaults

	util::options::Parser parser("< intervals", "> union_intervals");
	parser.parse(argv, argv + argc);

	try {
		IntervalList fIntervals, rIntervals, intervals;

		readIntervals(std::cin, fIntervals, rIntervals);
		
		unionIntervals(fIntervals);
		unionIntervals(rIntervals);

		mergeIntervals(fIntervals, rIntervals, intervals);

		writeIntervals(std::cout, intervals);

	} catch (const std::runtime_error& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
