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


#include "BasisBSplineMZ.hpp"
#include "BSpline.hpp"
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;


BasisBSplineMZ::
BasisBSplineMZ(vector<Basis*>& bases, DataMS& input, vector< vector<double> >& mzs, ii out_res, ii order, bool transient) :
	BasisBSpline(bases, 1, transient),
	is(input.get_is()),
	as(input.get_js().size()),
	ats(input.get_js().size())
{
	const vector<ii>& js = input.get_js();
	double mz_interval = 1.0033548378 / (60 * pow(2.0, (double) out_res));

	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	// find min and max m/z across spectra, and m for each A
	mz_min = DBL_MAX;
	mz_max = 0.0;
	vector<ii> ms(as.size());
	for (ii k = 0; k < (ii)as.size(); k++)
	{
		ms[k] = (ii)(is[k + 1] - is[k]);
		mz_min = mzs[js[k]].front() < mz_min ? mzs[js[k]].front() : mz_min;
		mz_max = mzs[js[k]].back() > mz_max ? mzs[js[k]].back() : mz_max;
	}

	// fill in b-spline mesh info
	mi.m = (ii) as.size();
	mi.ls[0] = out_res;
	mi.os[0] = (ii)floor(mz_min / mz_interval);
	mi.ns[0] = ((ii)ceil(mz_max / mz_interval)) + order - mi.os[0];

	// populate coo matrix
	ii done = 0;
	vector<ii> nnzs(as.size(), 0);
	BSpline bspline(order, 65536); // bspline basis function lookup table
	#pragma omp parallel for
	for (ii k = 0; k < (ii)as.size(); k++)
	{
		vector<fp> acoo;
		vector<ii> rowind;
		vector<ii> colind;

		for (ii i = 0; i < ms[k]; i++)
		{
			//if (gs[is[k] + i] >= 0.0)
			{
				double cf_min = mzs[js[k]][i] / mz_interval;
				double cf_max = mzs[js[k]][i + 1] / mz_interval;

				ii c_min = (ii)floor(cf_min);
				ii c_max = ((ii)ceil(cf_max)) + order;

				// work out basis coefficients
				for (ii c = c_min; c < c_max; c++)
				{
					double bf_min = (double)(c - order);
					double bf_max = (double)(c + 1);

					// intersection of bin and basis, between 0 and order+1
					double b_min = cf_min > bf_min ? cf_min - bf_min : 0.0;
					double b_max = cf_max < bf_max ? cf_max - bf_min : bf_max - bf_min;

					// basis coefficient b is _integral_ of area under b-spline basis
					fp b = (fp)(bspline.ibasis(b_max) - bspline.ibasis(b_min));

					if (b >= 0.000001) // small basis coefficients waste cpu and cause problems when calculating l2norm
					{
						acoo.push_back(b);
						rowind.push_back(i);
						colind.push_back(c - mi.os[0]);
					}
				}
			}
		}

		// create A
		as[k].init(ms[k], mi.ns[0], acoo, rowind, colind);
		ats[k].init(mi.ns[0], ms[k], acoo, colind, rowind);

		// display progress update
		#pragma omp critical
		{
			done++;
			if (done % 100 == 0)
			{
				for (int i = 0; i < 256; ++i) cout << '\b';
				cout << get_index() << " BasisBSplineMZ " << setw(1 + (int)(log10((float)mi.m))) << done << "/" << mi.m << " " << flush;
			}
		}

		nnzs[k] = (ii)acoo.size();
	}
	for (int i = 0; i < 256; ++i) cout << '\b';

	li nnz = 0; for (ii k = 0; k < (ii)as.size(); k++) nnz += nnzs[k];
	li mem = 0; for (ii k = 0; k < (ii)as.size(); k++) mem += as[k].mem();
	cout << get_index() << " BasisBSplineMZ";
	cout << " mz_range=(" << fixed << setprecision(2) << mz_min << ":" << setprecision(5) << mz_interval << ":" << setprecision(2) << mz_max << ")Th";
	cout << " C" << mi << " A[" << mi.n() << "," << is.back() << "]:" << nnz << " (" << (2 * mem) / 1024.0 / 1024.0 << "Mb)";
	if (transient) cout << " (t)";
	cout << endl;
}


void
BasisBSplineMZ::
synthesis(Matrix& f, const Matrix& c, bool accum) const
{
	for (ii k = 0; k < (ii) as.size(); k++)
	{
		//cout << "  synthesis BasisBsplineMZ C" << c << " x A" << ats[k] << " =";
		f.mult(ats[k], c, true, false);
		//cout << " F" << f << endl;
	}
}


void
BasisBSplineMZ::
analysis(Matrix& c_err, const Matrix& f_err, bool a_sqrd) const
{
	for (ii k = 0; k < (ii)as.size(); k++)
	{
		//cout << "  analysis BasisBsplineMZ Ferr" << f_err << " x A_t" << as[k] << " =";
		c_err.mult(as[k], f_err, false, a_sqrd);
		//cout << " Cerr" << c_err << endl;
	}
}

