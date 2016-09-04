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


#ifndef _SEAMASS_CORE_BASISBSPLINESCALE_HPP_
#define _SEAMASS_CORE_BASISBSPLINESCALE_HPP_


#include "BasisBSpline.hpp"


class BasisBSplineScale : public BasisBSpline
{ 
private:
	Matrix a;
	Matrix at;

	ii dim;

public:
	BasisBSplineScale(std::vector<Basis*>& bases, ii parent_index, ii dim, ii order = 3, bool transient = false);
	virtual ~BasisBSplineScale() {}

	void synthesis(Matrix& f, const Matrix& c, bool accum = true) const;
	void analysis(Matrix& c_err, const Matrix& f_err, bool a_sqrd = false) const;
};


#endif

