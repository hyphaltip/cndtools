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

#ifndef __BIO_ALIGNMENT_INTERVAL_HH__
#define __BIO_ALIGNMENT_INTERVAL_HH__

namespace bio { namespace alignment {

	struct Interval {
		size_t start;
		size_t end;

		Interval(const size_t start = 0, const size_t end = 0);
	};

	inline Interval::Interval(const size_t start, const size_t end)
		: start(start), end(end) {
	}
	
} }

#endif // __BIO_ALIGNMENT_INTERVAL_HH__
