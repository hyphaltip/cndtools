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

#ifndef __TYPES_HH__
#define __TYPES_HH__

#include <string>
using std::string;
#include <bitset>
using std::bitset;
#include <vector>
using std::vector;
#include <set>
using std::set;
#include <limits>
using std::numeric_limits;
#include <algorithm>
using std::sort;
using std::reverse;
using std::min;
using std::max;
using std::fill;
#include <functional>
using std::mem_fun;
using std::mem_fun_ref;
using std::bind2nd;
#include <utility>
using std::pair;
using std::make_pair;
#include <ostream>
using std::ostream;
#include <iterator>
using std::back_inserter;
#include <iostream>
using std::cout;
using std::cerr;
using std::clog;
using std::endl;
#include <cassert>

#include "filesystem.hh"
using namespace filesystem;

class Genome;
class Chromosome;
class Anchor;
class Edge;
class Clique;
class Run;

#define MAX_GENOMES 128

typedef bitset<MAX_GENOMES> Mask;
typedef long long GenomicDist;
typedef vector<Genome*>::const_iterator GenomeIter;

#endif // __TYPES_HH__
