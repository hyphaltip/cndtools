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

#include <cassert>

#include "bio/formats/maf/VariableMap.hh"

namespace bio { namespace formats { namespace maf {
	
	const std::string& VariableMap::getVariable(const std::string& var) const {
		assert(hasVariable(var));
		return variables.find(var)->second;
	}
	
	bool VariableMap::hasVariable(const std::string& var) const {
		return variables.find(var) != variables.end();
	}
	
	void VariableMap::setVariable(const std::string& var,
								  const std::string& val) {
		variables[var] = val;
	}

	VariableMap::~VariableMap() {}
	
} } }
