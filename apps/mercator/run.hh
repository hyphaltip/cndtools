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

#ifndef __RUN_HH__
#define __RUN_HH__

#include "types.hh"

class Run {
public:

	// Constructor
	Run();

	// Destructor - unclaims all cliques in this run (if they were
	// claimed by this run to begin with)
	~Run();

	bool hasAnchor(size_t num, size_t genome);

	// Set the minimum length (number of cliques) for a run to be
	// considered significant
	static void setMinRunLength(size_t num);

	// Get the minimum length (number of cliques) for a run to be
	// considered significant
	static size_t getMinRunLength();

	void setNum(const int n);
	int getNum() const;

	// Remove genomes from this run that only have one anchor.  This
	// is an iterative process in which genomes that originally had
	// more than one anchor may end up being removed if they were in
	// cliques with other removed genomes
	void removeSingletons();
	
	// Returns true if this run is significant
	bool isSignificant() const;

	// Set this run's visited status to VISITED
	void setVisited(bool visited = true);

	// Returns true if this run has already been visited during the
	// run joining procedure
	bool wasVisited() const;
	
	// Returns the mask indicating which genomes are contained in this
	// run
	Mask getMask() const;

	// Returns a vector containing the ordered cliques in this run
	const vector<Clique*>& getCliques() const;

	// Output the string representation of run R on the stream STRM
	friend ostream& operator<<(ostream& strm, const Run& r);

	// Output the string representation of this run on the standard
	// error stream
	void print() const;
	
	// Returns true if the anchors from the given GENOME are from
	// different chromosomes/anchors in this run.  This may only be the
	// case if this genome is a draft genome.
	bool isAssembled(const size_t genome) const;
	
	// Returns true if this run is to the left of the run OTHER,
	// adjacent to the other run in at least MINADJACENT genomes, and
	// within MAXDIST of OTHER in each genome shared between the two
	// runs.  If STRICT is set to true, the relevant ends of the two
	// runs that are on draft genomes need either be on the same
	// chromosome or be at the ends of their respective chromosomes.
	bool isLeftOf(const Run& other,
				  const size_t minAdjacent = 0,
				  const GenomicDist maxDist = numeric_limits<GenomicDist>::max(),
				  const bool strict = false) const;

	// Returns the next run from this run in GENOME in the direction
	// indicated by FORWARD.
	Run* nextRun(const size_t genome,
				 const bool forward = true) const;

	// Returns the run that is adjacent to the left side of this run
	// in GENOME
	Run* leftRun(const size_t genome) const;

	// Returns the run that is adjacent to the right side of this run
	// in GENOME
	Run* rightRun(const size_t genome) const;

	// Returns the last clique in this run
	Clique* getLastClique() const;

	// Returns the first clique in this run
	Clique* getFirstClique() const;

	// Returns in the vector GCLIQUES of all the cliques contained in
	// this run that are in the genome G
	void getCliques(vector<Clique*>& gcliques, size_t g) const;
	
	// Returns true if this run has no cliques
	bool isEmpty() const;

	// Add clique C to the end of this run
	void addCliqueToRight(Clique* c);

	// Append a run to the end of this run
	void addRunToRight(Run* r);

	// Split this run up into two runs by removing the right-most part
	// of this run that is on the same chromosome in GENOME.  Returns
	// the run that is split off the right.
	Run* splitRunRight(const size_t genome);

	// Split this run up into its component (one-clique) runs
	void split();
	
	// Returns the weight (score) of this run.  Currently this is the
	// sum of the weights of the cliques in this run.
	int getWeight() const;

	// Removes all cliques from this run.  Makes this run empty.
	void reset();

	// Flip the order and orientation of the anchors in this run.
	void flip();
	
	// Returns true if this run includes portions of the genomes
	// specified by GENOMEMASK
	bool hasGenomes(const Mask& genomeMask) const;

	// Returns true if this run includes portions of GENOME.
	bool hasGenome(const size_t genome) const;

	// Returns the number of cliques used to form this run
	size_t getLength() const;

	// Returns true if the anchors in the cliques in this run are in
	// increasing (coordinate) order in GENOME.  This function may
	// only be called if this run is not assembled for GENOME.
	bool isForward(const size_t genome) const;

	// Returns true if the right-most anchors in this run are in
	// increasing order in GENOME.  This is equivalent to isForward if
	// this run is not assembled for GENOME.
	bool isRightForward(const size_t genome) const;

	// Returns true if the left-most anchors in this run are in
	// increasing order in GENOME.  This is equivalent to isForward if
	// this run is not assembled for GENOME.
	bool isLeftForward(const size_t genome) const;

	// Returns...
	bool isRightChromEnd(const size_t genome) const;

	// Returns..
	bool isLeftChromEnd(const size_t genome) const;

	pair<Anchor*, Anchor*>
	getLeftInterRunAnchorInterval(const size_t genome) const;
	
	pair<Anchor*, Anchor*>
	getRightInterRunAnchorInterval(const size_t genome) const;

	void
	getLeftInterRunAnchorIntervals(vector< pair<Anchor*, Anchor*> >& ints) const;

	void
	getRightInterRunAnchorIntervals(vector< pair<Anchor*, Anchor*> >& ints) const;

	// Returns the left-most anchor in this run for GENOME.
	Anchor* getLeftAnchor(const size_t genome) const;

	// Returns the right-most anchor in this run for GENOME.
	Anchor* getRightAnchor(const size_t genome) const;

	// Returns the left-most clique in this run for GENOME.
	Clique* getLeftClique(const size_t genome) const;

	// Returns the right-most clique in this run for GENOME.
	Clique* getRightClique(const size_t genome) const;

	// Returns the chromosome of the left-most anchor in this run for
	// GENOME
	Chromosome* getLeftChrom(const size_t genome) const;

	// Returns the chromosome of the right-most anchor in this run for
	// GENOME
	Chromosome* getRightChrom(const size_t genome) const;
	
	// Returns the chromosome used in this run for GENOME.  May not be
	// called if this run is assembled for GENOME.
	Chromosome* getChrom(const size_t genome) const;

	GenomicDist leftCoord(const size_t genome) const;
	GenomicDist rightCoord(const size_t genome) const;
	GenomicDist leftSplitCoord(const size_t genome) const;
	GenomicDist rightSplitCoord(const size_t genome) const;
	GenomicDist startCoord(const size_t genome) const;
	GenomicDist endCoord(const size_t genome) const;	
	GenomicDist start(const size_t genome) const;
	GenomicDist end(const size_t genome) const;

	// Returns the number of cliques that are covered by this run in
	// genome GENOME (should this be the same as cliques.length?)
	size_t getNumCliques(const size_t genome) const;

	// Returns the number of anchors that are covered by this run in
	// genome GENOME
	size_t getNumAnchors(const size_t genome) const;

	// Returns the length of the region that is covered by this run in
	// genome GENOME
	GenomicDist getCoverage(const size_t genome) const;
	
	// Unclaims the cliques in this run, if they were claimed by this
	// run to begin with
	void unclaimCliques();

	// Claims the cliques in this run for itself and removes any other
	// runs that previously claimed the cliques in this run.  Also
	// sets the mask of the claimed cliques to be that of this run (?).
	void claimCliques();

	// Destroy all cliques inside of this run
	void deleteCliques();
	
	// Output a tab-delimited representation of this run (in terms of
	// the cliques that make it up) onto the stream STRM
	void dump(ostream& strm) const;

	// Output the intervals in each genome that are spanned by this
	// run. If EXTEND is true, the end points of this run will be
	// extended to the midpoints of the breakpoint regions on either
	// side of the run.
	void printMapLine(ostream& strm,
					  bool extend) const;

	// Output the pairwise hits between genomes within the cliques of
	// this run.
	void printPairwiseHits(ostream& strm) const;
	
	// Returns the run that is the result of joining this run to
	// OTHER, if possible, or NULL if it is not possible to join the
	// two runs.  OTHER must be to the left of this run.  MAXDIST is
	// the maximum distance allowed between adjacent anchors in the
	// run.
	Run* joinTo(Run& other, const GenomicDist maxDist);

	// Remove all edges incident to anchors internal to this run, but
	// that are not part of a clique
	void removeInterveningEdges();

	// Filter edges incident to anchors inside of runs
	void filterIntraRunEdges();

	// Filter edges incident to anchors to the left and right of this run
	void filterInterRunEdges();
	
private:
	static size_t minRunLength;
	
    Mask mask;
	vector<Clique*> cliques;
	vector<Clique*> starts;
	vector<Clique*> ends;

	// Flag variables for various algorithms acting on runs
	bool visited;
	bool joined;

	// Numbering for final output
	int num;
	
	// Returns true if this run has already been joined during the run
	// joining procedure
	bool wasJoined() const;

	// Returns true if this run can be joined to the run OTHER, to the
	// right of OTHER, within a maximum distance of MAXDIST in each
	// shared genome.  If it can be joined, the resulting joined run
	// is placed in JOINEDRUN.
	bool canJoinTo(Run& other, Run& joinedRun, const GenomicDist maxDist);

	// Reset the flags that get changed during the running of the
	// run joining procedure
	void resetJoinFlags();	
};

// Functor for comparing the lengths of runs
struct RunLengthSorter {
	bool operator()(const Run* r1, const Run* r2) const;
};

// Functor for comparing the coordinates of runs in GENOME
struct RunSorter {
	bool operator()(const Run* r1, const Run* r2) const;
};

// Order the runs in the vector RUNS according to their coordinates in
// GENOME
void orderRuns(vector<Run*>& runs);

// Return all of the runs in the vector RUNS
void getRuns(vector<Run*>& runs);

// Return all of the runs that contain at least one of the genomes in
// INCLUDED
void getRuns(vector<Run*>& runs, const Mask& included);

// Return all the runs that contain the genome G in the vector RUNS
void getRuns(vector<Run*>& runs, const size_t g);

// Join runs that are adjacent (on the same chromosome) in at least
// MINADJACENT genomes and can be joined without joining runs that are
// more than MAXDIST apart in any genome
void joinRuns(const size_t minAdjacent = 1,
			  const GenomicDist maxDist = numeric_limits<GenomicDist>::max());

// Number the runs in RUNS in order starting at 1
void numberRuns(vector<Run*>& runs);

void breakRuns();

void filterIntraRunEdges();

// Filter edges incident to anchors not inside of runs so that
// they are consistent with one of the two adjacent runs
void filterInterRunEdges();

void removeInsignificantRuns();
void removeSingletons();

// Break apart runs that form cycles in assembled genomes (i.e. if you
// traverse runs from left to right in an assembled genome and end up
// going in a circle)
void breakCyclicRuns();


// Implementation of inline Run methods
inline void Run::setNum(const int n) { num = n; }
inline int Run::getNum() const { return num; }

#endif // __RUN_HH__
