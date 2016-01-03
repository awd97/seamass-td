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


#ifndef _SEAMASS_CORE_CORE_HPP_
#define _SEAMASS_CORE_CORE_HPP_


#include <mkl.h>
#include <vector>


typedef float fp; // fp is the selected floating point precision (float or double)
typedef MKL_INT ii; // ii is the selected indexing integer size (32 or 64 bit)
typedef long long li; // 64bit integer


namespace utils
{
	void bin_mzs_intensities(std::vector< std::vector<double> >& mzs,
		                     std::vector< std::vector<double> >& intensities,
		                     ii instrument_type, std::vector<fp>& exposures);

	void create_gs(std::vector<fp>& gs, std::vector<li>& is, std::vector<ii>& js,
		           const std::vector< std::vector<double> >& intensities);
}


#endif

