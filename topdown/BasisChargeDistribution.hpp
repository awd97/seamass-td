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


#ifndef _SEAMASS_TOPDOWN_BASISCHARGEDISTRIBUTION_HPP_
#define _SEAMASS_TOPDOWN_BASISCHARGEDISTRIBUTION_HPP_


#include "../core/Basis.hpp"
#include "../core/BSpline.hpp"


class BasisChargeDistribution : public Basis
{
protected:
	// Input (gs) data
	const std::vector<li>& is;

	// Sparse basis matrices A (one per input spectrum)
	std::vector<SparseMatrix> as;

	// Output (cs) metadata
	li nc; // number of coefficients across all spectra
	std::vector<ii> ci0s, ci1s, cos; // coefficient range and offset for each charge state 
	double mass_interval; // interval between output coefficients in Daltons

	// Helpers
	static double proton_mass; // mass of a proton (positive charge) in Daltons
	static double peak_width(double fwhm, double mz); // orbitrap peak width given fwhm @ 400 m/z 

public:
	BasisChargeDistribution(std::vector<Basis*>& bases,
		const std::vector< std::vector<double> >& mzs,
		const std::vector<fp>& gs, const std::vector<li>& is, const std::vector<ii>& js,
		double out_interval, ii out_res, ii factor_res, ii max_z, double peak_fwhm,
		bool transient = false);

	~BasisChargeDistribution();

	void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true) const;
	void analysis(std::vector<fp>& es, const std::vector<fp>& fs) const;
	void l2norm(std::vector<fp>& es, const std::vector<fp>& fs) const;
	void shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage) const;
	li get_nc() const { return nc; }

	const std::vector<ii>& get_ci0s() const { return ci0s; }
	const std::vector<ii>& get_ci1s() const { return ci1s; }
	const std::vector<ii>& get_cos() const { return cos; }
	ii get_ns() const { return as.size(); }
	double get_mass_interval() const { return mass_interval; }

	void write_cs(const std::vector<fp>& cs) const;
};


#endif

