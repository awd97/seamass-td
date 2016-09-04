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


#ifndef _SEAMASS_CORE_DATAMS_HPP_
#define _SEAMASS_CORE_DATAMS_HPP_


#include "Matrix.hpp"


class DataMS
{
private:
	std::vector<fp> exposures; // per-spectrum exposures
	Matrix g;          // raw data as a 1D vector
	std::vector<ii> is; // spectrum indicies into gs, size is number of spectra + 1 
	std::vector<ii> js; // spectrum indicies into mzs & intensitives, size is number of spectra

public:
	DataMS(const std::string& id, const std::string& config_id, int instrument_type,
		std::vector<double>& sts, std::vector< std::vector<double> >& mzs, std::vector< std::vector<double> >& intensities);
	~DataMS() {}

	const Matrix& get_g() const { return g; }
	const std::vector<ii>& get_is() const { return is; }
	const std::vector<ii>& get_js() const { return js; }
};


#endif

