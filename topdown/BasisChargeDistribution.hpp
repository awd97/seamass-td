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


#ifndef _SEAMASS_BASISCHARGEDISTRIBUTION_HPP_
#define _SEAMASS_BASISCHARGEDISTRIBUTION_HPP_


#include "../core/Basis.hpp"
#include "../core/BSpline.hpp"


class BasisUniformChargeDistribution : public Basis
{
protected:
    const std::vector<li>& is;
    
    // CSR sparse A basis matrices and their transposes
    std::vector<ii> nnz;
    std::vector<ii> m;
    std::vector< std::vector<fp> > a;
    std::vector< std::vector<ii> > ia, ja;
 	
	double rc; // resolution of output (in Daltons)
    
public:
	BasisUniformChargeDistribution(std::vector<Basis*>& bases,
                            const std::vector< std::vector<double> >& mzs,
                            const std::vector<fp>& gs, const std::vector<li>& is, const std::vector<ii>& js,
							ii mass_res, ii max_z, double max_peak_width,
							bool transient = false);
    
	~BasisUniformChargeDistribution();
    
	void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true);
	void analysis(std::vector<fp>& es, const std::vector<fp>& fs);
	void l2norm(std::vector<fp>& es, const std::vector<fp>& fs);

	void write_cs(const std::vector<fp>& cs);
};


class BasisFreeformChargeDistribution : public Basis
{
protected:
	const std::vector<li>& is;

	// CSR sparse A basis matrices and their transposes
	std::vector<ii> nnz;
	std::vector<ii> m;
	std::vector< std::vector<fp> > a;
	std::vector< std::vector<ii> > ia, ja;

	double rc; // resolution of output (in Daltons)
	ii n, max_z;

public:
	BasisFreeformChargeDistribution(std::vector<Basis*>& bases,
		const std::vector< std::vector<double> >& mzs,
		const std::vector<fp>& gs, const std::vector<li>& is, const std::vector<ii>& js,
		ii mass_res, ii max_z, double max_peak_width,
		bool transient = false);

	~BasisFreeformChargeDistribution();

	void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true);
	void analysis(std::vector<fp>& es, const std::vector<fp>& fs);
	void l2norm(std::vector<fp>& es, const std::vector<fp>& fs);
	void shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage);

	void write_cs(const std::vector<fp>& cs);
};


#endif

