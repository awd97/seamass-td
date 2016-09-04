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


#ifndef _SEAMASS_CORE_BASISBSPLINE_HPP_
#define _SEAMASS_CORE_BASISBSPLINE_HPP_


#include "Basis.hpp"


class BasisBSpline : public Basis
{
public:
	struct MeshInfo
	{
		ii d;				 // dimension of each b-spline coefficients mesh
		std::vector<ii> ls;  // dyadic level for each dimension
		std::vector<ii> os;  // coefficient offset for each dimension
		std::vector<ii> ns;  // number of coefficients for each dimension (make up the columns)
		ii m;                // number of meshes (columns)

		MeshInfo(ii d);
		~MeshInfo();

		void operator=(const MeshInfo& bo);
		ii n() const;
		li size() const;
	};

protected:
	MeshInfo mi;

public:
	BasisBSpline(std::vector<Basis*>& bases, ii dims, bool transient = false, ii parent_index = -1);
	virtual ~BasisBSpline() {}

	li get_nc() const { return mi.size(); }

	const MeshInfo get_meshinfo() const { return mi; }
};

std::ostream& operator<<(std::ostream& os, const BasisBSpline::MeshInfo& mi);

#endif

