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


#include "BasisChargeDistribution.hpp"
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;


double
BasisChargeDistribution::
proton_mass = 1.007276466879; // in Daltons

double
BasisChargeDistribution::
peak_width(double fwhm, double mz)
{
	return (400.0 / (sqrt(400.0 / mz) * fwhm)) / 0.361175; // fwhm of cubic basis function of unit width is 0.361175 
}


BasisChargeDistribution::
BasisChargeDistribution(vector<Basis*>& bases,
                        const vector< vector<double> >& mzs,
						const vector<fp>& gs, const vector<li>& _is, const vector<ii>& js,
						double _mass_interval, ii out_res, ii factor_res, ii max_z, double peak_fwhm,
						bool transient) :
	Basis(bases, 0, transient),
	is(_is),
	as(js.size()),
	ci0s(max_z, numeric_limits<ii>::max()),
	ci1s(max_z, 0),
	cos(max_z+1),
	mass_interval(_mass_interval / pow(2.0, (double)out_res))
{
	cout << get_index() << " BasisChargeDistribution " << flush;

	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix
	vector<ii> ms(js.size());
	vector<ii> nnzs(js.size(), 0);

	// find min (mz0) and max (mz1) m/z across all spectra
	double mz0 = numeric_limits<double>::max();
	double mz1 = 0.0;
	for (ii j = 0; j < js.size(); j++)
	{
		ms[j] = (ii)(is[j + 1] - is[j]);
		mz0 = mzs[js[j]].front() < mz0 ? mzs[js[j]].front() : mz0;
		mz1 = mzs[js[j]].back() > mz1 ? mzs[js[j]].back() : mz1;
	}
	mz0 -= 0.5*peak_width(peak_fwhm, mz0);
	mz1 += 0.5*peak_width(peak_fwhm, mz1);

	// calculate output mass and hence coefficient range
	double mass0 = mz0 - proton_mass;
	double mass1 = max_z * (mz1 - proton_mass);
	ii mi0 = ceil(mass0 / mass_interval);
	ii mi1 = floor(mass1 / mass_interval);
	mass0 = mi0 * mass_interval;
	mass1 = mi1 * mass_interval;

	// figure out nnz, ci0s and ci1s
	#pragma omp parallel for
	for (ii j = 0; j < js.size(); ++j)
	{
		for (ii z = 0; z < max_z; z++)
		for (ii i = 0; i < ms[j]; i++)
		if (gs[is[j] + i] >= 0.0)
		{
			double cf0 = (z + 1) * (mzs[js[j]][i] - 0.5*peak_width(peak_fwhm, mzs[js[j]][i]) - proton_mass) / mass_interval;
			double cf1 = (z + 1) * (mzs[js[j]][i + 1] + 0.5*peak_width(peak_fwhm, mzs[js[j]][i + 1]) - proton_mass) / mass_interval;

			ii ci0 = (ii)ceil(cf0);
			ii ci1 = (ii)floor(cf1);

			nnzs[j] += ci1 - ci0 + 1;
			ci0s[z] = ci0s[z] < ci0 ? ci0s[z] : ci0;
			ci1s[z] = ci1s[z] > ci1 ? ci1s[z] : ci1;
		}
	}

	li nc_per_spectrum = 0;
	cos[0] = 0;
	for (ii z = 0; z < max_z; z++)
	{
		cos[z+1] = cos[z] + ci1s[z] - ci0s[z] + 1;
	}
	nc = cos[max_z] * js.size();

	// populate coo matrix
	ii done = 0;
	BSpline bspline(3, 65535); // bspline basis function lookup table
	#pragma omp parallel for
	for (ii j = 0; j < js.size(); ++j)
	{
		vector<fp> acoo(nnzs[j]);
		vector<ii> rowind(nnzs[j]);
		vector<ii> colind(nnzs[j]);

		ii k = 0;
		for (ii z = 0; z < max_z; z++)
		{
			for (ii i = 0; i < ms[j]; i++)
			{
				if (gs[is[j] + i] >= 0.0)
				{
					// compute the range of peak coefficient indicies [ci0, ci1] that overlap with the raw data bin [mzs[js[j]][i], mzs[js[j]][i+1]]   
					double cf0 = (z + 1) * (mzs[js[j]][i] - 0.5*peak_width(peak_fwhm, mzs[js[j]][i]) - proton_mass) / mass_interval;
					double cf1 = (z + 1) * (mzs[js[j]][i + 1] + 0.5*peak_width(peak_fwhm, mzs[js[j]][i + 1]) - proton_mass) / mass_interval;
					ii ci0 = (ii)ceil(cf0);
					ii ci1 = (ii)floor(cf1);

					// for each coefficient index, determine what proportion of that peak overlaps with the raw data bin
					// (Gaussian peak shape is approximated by a cubic b-spline basis function) 
					for (ii ci = ci0; ci <= ci1; ci++)
					{
						double basisM = (ci * mass_interval) / (z + 1) + proton_mass;
						double pw = peak_width(peak_fwhm, basisM);
						double basis0 = basisM - 0.5*pw;
						double basis1 = basisM + 0.5*pw;

						// intersection of bin and basis
						double b0 = mzs[js[j]][i] > basis0 ? (mzs[js[j]][i] - basis0) / (basis1 - basis0) : 0.0;
						double b1 = mzs[js[j]][i + 1] < basis1 ? (mzs[js[j]][i + 1] - basis0) / (basis1 - basis0) : 1.0;

						// basis coefficient b is _integral_ of area under b-spline basis
						double b = bspline.ibasis(b1) - bspline.ibasis(b0);
						if (b <= FLT_MIN) b = FLT_MIN;

						acoo[k] = b;
						rowind[k] = i;
						colind[k] = cos[z] + (ci - ci0s[z]);

						k++;
					}
				}
			}

			// display progress update
			#pragma omp critical
			{
				done++;
				if (done % 1 == 0)
				{
					for (int i = 0; i < 256; ++i) cout << '\b';
					cout << get_index() << " BasisChargeDistribution " << setw(1 + (int)(log10((float)js.size()*max_z))) << done << "/" << js.size()*max_z << " " << flush;
				}
			}
		}

		// create As
		as[j].init(ms[j], cos[max_z], acoo, rowind, colind);
	}
	for (int i = 0; i < 256; ++i) cout << '\b';

	cout << get_index() << " BasisChargeDistribution";
	cout << " mass_range=[" << setprecision(2) << mass0 << ":" << setprecision(4) << mass_interval << ":" << setprecision(2) << mass1 << "]Da";
	if (transient) cout << " (t)" << endl; else cout << " nc=" << nc << endl;
	for (ii j = 0; j < js.size(); j++)
	{
		cout << endl;
		cout << "  A[" << j << "]=";
		as[j].print(cout);
	}

	for (ii z = 0; z < max_z; z++)
	{
		cout << endl << "  z=" << z+1 << " co=" << cos[z] << " ci_range=[" << ci0s[z] << ":" << ci1s[z] << "]" << flush;
	}
	cout << endl;
}


BasisChargeDistribution::~BasisChargeDistribution()
{
}


void
BasisChargeDistribution::
synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum) const
{
	for (li j = 0; j < as.size(); j++)
	{
		as[j].mult(&(fs.data()[is[j]]), &(cs.data()[j*as[j].get_n()]), false, accum);
	}
}


void
BasisChargeDistribution::
analysis(vector<fp>& es, const vector<fp>& fs) const
{
	for (li j = 0; j < as.size(); j++)
	{
		as[j].mult(&(es.data()[j*as[j].get_n()]), &(fs.data()[is[j]]), true);
	}
}


void
BasisChargeDistribution::
l2norm(vector<fp>& es, const vector<fp>& fs) const
{
	for (li j = 0; j < as.size(); j++)
	{
		as[j].sqr_mult(&(es.data()[j*as[j].get_n()]), &(fs.data()[is[j]]), true);
	}
}


void
BasisChargeDistribution::
shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage) const
{
	// GROUP-WISE SHRINKAGE! (note - intuitive implemention, not mathematically verified yet)
	// also horrible gather operation
	for (li j = 0; j < as.size(); j++)
	{
		vector<fp> gcs(ci1s.back() - ci0s.front() + 1, 0.0);

		// sum up coefficients per group
		for (ii z = 0; z < cos.size() - 1; z++)
		for (ii i = cos[z]; i < cos[z + 1]; i++)
		if (es[j*cos.back() + i] > 0.0)
		{
			gcs[i - cos[z] + ci0s[z] - ci0s.front()] += cs[j*cos.back() + i];
		}

		// scale the shrinkage to be proportional to the contribution of this coefficient to the group total
		#pragma omp parallel for
		for (ii z = 0; z < cos.size(); z++)
		for (ii i = cos[z]; i < cos[z + 1]; i++)
		if (es[j*cos.back() + i] > 0.0)
		{
			if (cs[j*cos.back() + i] > 0.0)
			{
				double scale = cs[j*cos.back() + i] / gcs[i - cos[z] + ci0s[z] - ci0s.front()];
				es[j*cos.back() + i] *= cs[j*cos.back() + i] / (scale * shrinkage * l2[j*cos.back() + i] + wcs[j*cos.back() + i]);
			}
			else
			{
				es[j*cos.back() + i] = 0.0;
			}
		}
	}
}


void
BasisChargeDistribution::
write_cs(const std::vector<fp>& cs) const
{
	for (li j = 0; j < as.size(); j++)
	{
		ostringstream oss; oss << "profile" << j << ".csv";
		ofstream ofs(oss.str());

		// sum up coefficients per group
		vector<fp> gcs(ci1s.back() - ci0s.front() + 1, 0.0);
		for (ii z = 0; z < cos.size() - 1; z++)
		for (ii i = cos[z]; i < cos[z + 1]; i++)
		{
			gcs[i - cos[z] + ci0s[z] - ci0s.front()] += cs[j*cos.back() + i];
		}

		ofs << "mass,intensity" << setprecision(10) << endl;
		for (ii i = 0; i < gcs.size(); ++i)
		{
			if (gcs[i] > 0.0 || i > 0 && gcs[i - 1] > 0.0 || i < gcs.size() - 1 && gcs[i + 1] > 0.0)
			{
				ofs << (ci0s.front() + i) * mass_interval << "," << gcs[i] << endl;
			}
		}
	}
}
