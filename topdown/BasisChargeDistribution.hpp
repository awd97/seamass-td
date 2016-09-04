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


#include "../core/BasisBSpline.hpp"


/*class BasisChargeDistribution : public Basis
{
protected:
	// Input (gs) data indicies
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

	void synthesis(SparseMatrix& fs, const SparseMatrix& cs, bool accum = true) const;
	void analysis(SparseMatrix& es, const SparseMatrix& fs) const;
	void l2norm(SparseMatrix& es, const SparseMatrix& fs) const;
	void shrink(SparseMatrix& es, const SparseMatrix& cs, const SparseMatrix& l2, const SparseMatrix& wcs, double shrinkage) const;
	li get_nc() const { return nc; }

	const std::vector<ii>& get_ci0s() const { return ci0s; }
	const std::vector<ii>& get_ci1s() const { return ci1s; }
	const std::vector<ii>& get_cos() const { return cos; }
	ii get_ns() const { return (ii) as.size(); }
	double get_mass_interval() const { return mass_interval; }

	void write_cs(const std::vector<fp>& cs) const;
};*/


class BasisBSplineMZ : public BasisBSpline
{
protected:
	// Input (gs) data indicies
	const std::vector<li>& is;

	// Sparse basis matrices A (one per input spectrum)
	std::vector<SparseMatrix> as;

	// min and max m/z across all spectra
	double mz_min, mz_max;

public:
	BasisBSplineMZ(std::vector<Basis*>& bases, const std::vector< std::vector<double> >& mzs,
		const std::vector<fp>& gs, const std::vector<li>& is, const std::vector<ii>& js,
		ii out_res, ii order = 3,
		bool transient = false);

	~BasisBSplineMZ();

	void synthesis(SparseMatrix& fs, const SparseMatrix& cs, bool accum = true) const;
	void analysis(SparseMatrix& es, const SparseMatrix& fs) const;
	void l2norm(SparseMatrix& es, const SparseMatrix& fs) const;

	double get_min() const { return mz_min; }
	double get_max() const { return mz_max; }
};


#endif

