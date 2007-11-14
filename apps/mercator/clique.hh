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

#ifndef __CLIQUE_HH__
#define __CLIQUE_HH__

#include "types.hh"

class Clique {
public:
	// Make a new clique
	Clique();

	// Destroy this clique.  Any anchors in this clique are unclaimed.
	~Clique();

	// Mark this clique as one to keep in the final map
	void keepClique();

	// Returns true if this clique has been marked as one to keep
	bool isKept() const;
	
	// Flip the orientations (strands) of all the anchors in this clique
	void flip();
	
	// Sets the genomes that are members of this clique to be those
	// specified by NEWMASK.  If this clique originally had anchors
	// from genomes that are not part of NEWMASK, those anchors are
	// unclaimed and removed from this clique.
	void setMask(const Mask& newMask);

	// Makes this clique a member of the run R
	void setRun(Run* r);

	// Add the anchor A to this clique
	void addAnchor(Anchor* a);

	// Remove anchor A from this clique
	void removeAnchor(Anchor* a);

	// remove anchor from genome G from this clique
	void removeAnchor(size_t g);

	// Unclaims the anchors in this clique
	void removeAnchorPtrs();

	// Claims the anchors in this clique and removes all edges
	// incident to these anchors that are not part of the clique
	void useClique();

	// Returns the mask for this clique that specifies which genomes
	// have anchors in this clique
	Mask getMask() const;

	// Returns the anchor in this clique for genome num G
	Anchor* getAnchor(size_t g) const;

	// Returns the anchor in this clique for genome G
	Anchor* getAnchor(Genome* g) const;

	// Returns the run that this clique is a member of or NULL if this
	// clique is not a member of any run
	Run* getRun() const;

	// Returns the number of genomes that have anchors in this clique
	size_t getSize() const;

	// Returns the weight (score) of this clique
	int getWeight(const Mask& weightMask) const;

	// Returns the next clique in the genome G in the direction
	// indicated by FORWARD
	Clique* nextClique(size_t g, bool forward = true);

	// Returns the next run in the genome G in the direction indicated
	// by FORWARD
	Run* nextRun(size_t g, bool forward = true);
	
	// Returns true if this is the first clique on the chromosome of
	// the genome G
	bool isFirstClique(const size_t g) const;

	// Returns true if this is the last clique on the chromosome of
	// the genome G
	bool isLastClique(const size_t g) const;

	// Returns true if this is the first or last clique on the
	// chromosome of genome G
	bool isEndClique(const size_t g) const;

	// Returns true if this clique is a member of a run
	bool isInRun() const;

	// Returns true if this clique includes an anchor from genome G
	bool hasGenome(size_t g) const;

	// Returns true if this clique includes anchors from the genomes
	// specified by GENOMEMASK
	bool hasGenomes(const Mask& genomeMask) const;
	
	// Removes all cliques in between this clique and OTHER in the
	// genomes specified by REMOVEMASK
	void removeCliquesInBetween(const Clique* other,
								const Mask removeMask);


	// Removes all edges incident to anchors in between this clique
	// and OTHER in GENOME
	void removeEdgesInBetween(const Clique* other, const size_t genome);
	
	// Output stream operator for debugging cliques
	friend ostream& operator<<(ostream& strm, const Clique& c);

	bool isConnected() const;
	
	void print() const;

	friend bool operator==(const Clique& c1, const Clique& c2);
	
private:
    Mask mask;
	vector<Anchor*> anchors;
	Run* run;
	bool keep;
};

#endif // __CLIQUE_HH__
