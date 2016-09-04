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


#include "BasisBSplineScale.hpp"
#include "BSpline.hpp"
#include <limits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;



BasisBSplineScale::
BasisBSplineScale(vector<Basis*>& bases, ii parent_index, ii _dim, ii order, bool transient) :
BasisBSpline(bases, static_cast<BasisBSpline*>(bases[parent_index])->get_meshinfo().d, transient, parent_index),
	dim(_dim)
{
	const MeshInfo parent_mi = static_cast<BasisBSpline*>(bases[parent_index])->get_meshinfo();
	mi = parent_mi;
	mi.ls[dim] = parent_mi.ls[dim] - 1;
	mi.os[dim] = parent_mi.os[dim] / 2;
	mi.ns[dim] = (parent_mi.os[dim] + parent_mi.ns[dim] - 1 - order) / 2 + order + 1 - mi.os[dim];
	ii m = parent_mi.ns[dim];
	ii n = mi.ns[dim];

	ii stride = 1;
	for (ii j = 0; j < dim; j++) stride *= mi.ns[j];

	// create our kernel
	ii nh = order + 2;
	vector<fp> hs(nh);
	double sum = 0.0;
	for (ii i = 0; i < nh; i++)
	{
		hs[i] = (fp) (1.0 / pow(2.0, (double)order) * BSpline::factorial(order + 1) / (double)(BSpline::factorial(i)*BSpline::factorial(order + 1 - i)));
		sum += hs[i];
	}
	for (ii i = 0; i < nh; i++)
	{
		hs[i] /= (fp) sum;
	}

	// create A as a temporary COO matrix
	vector<fp> acoo(nh * n);
	vector<ii> rowind(nh * n);
	vector<ii> colind(nh * n);

	ii nnz = 0;
	ii offset = order + ((parent_mi.os[dim] + 1) % 2);
	for (ii j = 0; j < n; j++)
	{
		for (ii i = 0; i < nh; i++)
		{
			rowind[nnz] = 2 * j + i - offset;
			if (rowind[nnz] < 0 || rowind[nnz] >= m) continue;
			acoo[nnz] = hs[i];
			colind[nnz] = j;

			nnz++;
		}
	}

	a.init(m, n, acoo, rowind, colind);
	at.init(n, m, acoo, colind, rowind);

	cout << get_index() << " BasisBSplineScale parent=" << get_parent_index() << " dim=" << dim << " C" << mi;
	cout << " A" << a << " (" << (a.mem() + at.mem()) / 1024.0 / 1024.0 << "Mb)";
	if (transient) cout << " (t)";
	cout << endl;
}


void
BasisBSplineScale::
synthesis(Matrix& f, const Matrix& c, bool accum) const
{
	//cout << "  synthesis BasisBsplineScale C" << c << " x A" << at << " =";
	f.mult(at, c, accum, false);
	//cout << " F" << f << endl;
}


void
BasisBSplineScale::
analysis(Matrix& c_err, const Matrix& f_err, bool a_sqrd) const
{
	//cout << "  analysis BasisBSplineScale Ferr" << f_err << " x A_t" << a << " =";
	c_err.mult(a, f_err, false, a_sqrd);
	//cout << " Cerr" << c_err << endl;
}
