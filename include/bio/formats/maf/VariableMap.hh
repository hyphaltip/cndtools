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

#ifndef __BIO_FORMATS_MAF_VARIABLEMAP_HH__
#define __BIO_FORMATS_MAF_VARIABLEMAP_HH__

#include <string>

#include "util/stl.hh"

namespace bio { namespace formats { namespace maf {

    struct VariableMap {
		typedef util::stl::hash_map<std::string, std::string> Map;
		Map variables;

		const std::string& getVariable(const std::string& var) const;
		bool hasVariable(const std::string& var) const;
		void setVariable(const std::string& var,
						 const std::string& val);

		virtual ~VariableMap();
	};

} } }

#endif // __BIO_FORMATS_MAF_VARIABLEMAP_HH__
