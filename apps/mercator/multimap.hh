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

#ifndef __MULTIMAP_HH__
#define __MULTIMAP_HH__

#include "types.hh"

struct MapOptions {
	int repeatNum;
	float repeatPct;
	double maxE;
	float prunePct;
	int maxDist;
	int minRunLength;
	int padding;
	std::string indir;
	std::string outdir;
	bool quiet;
	bool outputHits;
	bool outputRuns;
	std::string phitFilename;
};

void removeAllCliques();

void removeUnusedCliques();

void removeRuns();

int countCliquesInRuns(vector<Run*>& runs);

int countCliquesInGenome(size_t genome);

void printNonUsedCliques(ostream& strm,
						 size_t genome);

void findRunsForMask(const Mask& mask,
					 const GenomicDist maxDist,
					 bool debug = false);
	

void findRuns(const vector<Mask>& masks,
			  const GenomicDist maxDist,
			  vector<Run*>& runs,
			  bool debug = false);

void printRuns(const vector<Run*>& runs,
			  ostream& strm);

void printMap(const vector<Run*>& runs,
			  ostream& strm,
			  const bool extend);

void printPairwiseHits(const vector<Run*>& runs,
					   ostream& strm);

size_t findCliques(const Mask m);

void findCliques(const vector<Mask>& masks);

struct RunStats {
	int numRuns;
	int numCliques;
	double meanRunLength;
	double medianRunLength;
	double maxRunLength;
	double minRunLength;
};

RunStats calcRunStats(vector<Run*>& runs);

void makeMap(vector<Run*>& runs, const MapOptions& options);
void joinPairwiseMaps(vector<Run*>& runs, const MapOptions& options);
void joinPairwiseMaps2(vector<Run*>& runs, const MapOptions& options);

#endif // __MULTIMAP_HH__
