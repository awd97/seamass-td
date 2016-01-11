//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  biospi Laboratory, EEE, University of Liverpool, UK
//
// This file is part of seaMass-TD.
//
// seaMass-TD is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// seaMass is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with seaMass-TD.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef _SEAMASS_TOPDOWN_TOPDOWN_HPP_
#define _SEAMASS_TOPDOWN_TOPDOWN_HPP_


//#include "seamass_topdown_export.h"
#include <vector>
#include <string>


namespace seamass
{
	//SEAMASS_TOPDOWN_EXPORT
	void topdown_notice();

	//SEAMASS_TOPDOWN_EXPORT
	void topdown(const std::string& id,
				 const std::string& config_id,
				 int instrument_type,
				 std::vector<double>& rts,
				 std::vector< std::vector<double> >& mzs,
				 std::vector< std::vector<double> >& intensities,
				 double out_mass0, double out_mass1, int out_res, int max_z, double max_peak_width, int shrink, int tol,
				 int threads = 1, int debug = 0);
}


#endif
