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

#ifndef __POLYALIGN_HH__
#define __POLYALIGN_HH__

#include <string>
#include <vector>
#include <map>

#include "polytope/Polytope.hh"
#include "bio/alignment/DNAScoringMatrix.hh"

typedef polytope::Polytope<int> IntegerPolytope;
typedef std::vector<IntegerPolytope> PolytopeList;
typedef bio::alignment::DNAScoringMatrix<IntegerPolytope> PolytopeMatrix;

typedef std::string Variable;
typedef std::set<Variable> VariableSet;
typedef std::map<Variable, IntegerPolytope> VariablePolytopeMap;
typedef std::map<Variable, size_t> VariableDimensionMap;
typedef bio::alignment::DNAScoringMatrix<Variable> VariableMatrix;

typedef bio::alignment::DNAScoringMatrix<int> IntegerMatrix;

IntegerPolytope calculatePolytope(const std::string& seq1,
								  const std::string& seq2,
								  PolytopeMatrix& matrix,
								  IntegerPolytope& space,
								  IntegerPolytope& gap,
								  size_t numVars);

void calcMiddleRowPolytopes(const std::string& seq1,
							const std::string& seq2,
							size_t num_params,
							PolytopeList& forwardMatch,
							PolytopeList& forwardGap1,
							PolytopeList& forwardGap2,
							PolytopeList& backwardMatch,
							PolytopeList& backwardGap1,
							PolytopeList& backwardGap2);

extern const std::string ZERO;

VariableSet
makeVariableSet(VariableMatrix& matrix, Variable& space, Variable& gap);

size_t
getNumVars(const VariableSet& varSet);

VariablePolytopeMap
makeVariablePolytopeMap(VariableSet& varSet);

PolytopeMatrix
makePolytopeMatrix(VariableMatrix& varMatrix,
				   VariablePolytopeMap& varMap);

VariableDimensionMap
makeVariableDimensionMap(VariableSet& varSet);

void
setParameters(IntegerMatrix& scoreMatrix, int& spaceScore, int& gapScore,
			  VariableMatrix& varMatrix, Variable& space, Variable& gap,
			  VariableSet& varSet,
			  IntegerPolytope::Ray& paramRay);

#endif // __POLYALIGN_HH__
