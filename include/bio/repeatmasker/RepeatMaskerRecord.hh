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

#ifndef __BIO_REPEATMASKER_REPEATMASKERRECORD_HH__
#define __BIO_REPEATMASKER_REPEATMASKERRECORD_HH__

#include <vector>

#include "bio/genome/BasicInterval.hh"

namespace bio { namespace repeatmasker {

    class RepeatMaskerRecord {
	public:
		RepeatMaskerRecord();

		int getScore() const;
		float getPctDivergence() const;
		float getPctDeleted() const;
		float getPctInserted() const;
		std::string getQueryName() const;
		genome::Position getQueryStart() const;
		genome::Position getQueryEnd() const;
		genome::BasicInterval getQueryInterval() const;
		genome::Distance getQueryLeft() const;
		bool isComplement() const;
		std::string getRepeatClass() const;
		std::string getRepeatName() const;
		genome::Position getRepeatStart() const;
		genome::Position getRepeatEnd() const;
		genome::BasicInterval getRepeatInterval() const;
		genome::Distance getRepeatLeft() const;
		int getID() const;
		bool isIncludedInHigherScoringMatch() const;

		friend std::ostream& operator<<(std::ostream& strm,
										const RepeatMaskerRecord& r);

		friend std::istream& operator>>(std::istream& strm,
										RepeatMaskerRecord& r);
		
		friend void operator>>(const std::string& line,
							   RepeatMaskerRecord& r);
		
	private:
		std::vector<std::string> fields;
    };

} }

#endif // __BIO_REPEATMASKER_REPEATMASKERRECORD_HH__
