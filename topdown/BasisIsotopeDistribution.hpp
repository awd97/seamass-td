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


#ifndef _SEAMASS_BASISISOTOPEDISTRIBUTION_HPP_
#define _SEAMASS_BASISISOTOPEDISTRIBUTION_HPP_


#include "../core/Basis.hpp"
#include "../core/BSpline.hpp"


class BasisIsotopeDistribution : public Basis
{
protected:
	// CSR sparse A basis matrix
	ii m, nnz;
	std::vector<fp> a;
	std::vector<ii> ia, ja;

public:
	BasisIsotopeDistribution(std::vector<Basis*>& bases, Basis* parent,
		                     ii out_res, ii factor_res, ii max_z,
							 bool transient = false);

	~BasisIsotopeDistribution();

	void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true);
	void analysis(std::vector<fp>& es, const std::vector<fp>& fs);
	void l2norm(std::vector<fp>& es, const std::vector<fp>& fs);
	void shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage);
};


#endif

