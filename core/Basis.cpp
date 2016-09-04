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


#include "Basis.hpp"
#include <limits>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;


Basis::
Basis(vector<Basis*>& bases, bool _transient, ii _parent_index) :
	parent_index(_parent_index),
	transient(_transient)
{
	index = (ii) bases.size();
	bases.push_back(this);
	if (parent_index >= 0) bases[parent_index]->child_count++;
}


Basis::ErrorInfo
Basis::
error(Matrix& f_err, const Matrix& f, const Matrix& g) const
{
	//cout << "  error Basis G" << g << " %/% F" << f << " =";
	f_err.elem_div(f, g);
	//cout << " F_err" << f_err << endl;

	/*if (fs[i] > 0.0 && gs[i] >= 0.0)
	{
		fs[i] = gs[i] / fs[i];
	}
	else
	{
		fs[i] = 1.0;
	}*/


	/*ErrorInfo info;

	/* ii size_d = 0;
	double dis = 0.0;
	double err = 0.0;
	double sum_g = 0.0;
	double sum_f = 0.0;
	info.max_error = 0;
	for (li i = 0; i < gs.size(); i++)
	{
	double v = fabs(gs[i]-fs[i]);
	info.max_error = info.max_error > v ? info.max_error : v;
	}
	// bug in this openmp section at present for err and discrep
	//#pragma omp parallel for simd reduction(+:dis,err,size_d,sum_g,sum_f)
	for (li i = 0; i < gs.size(); i++)
	{
	sum_g += gs[i];
	sum_f += fs[i];

	if (fs[i] > 0.0 && gs[i] >= 0.0)
	{
	dis += ((gs[i]-ceil(fs[i]))*(gs[i]-ceil(fs[i])))/ceil(fs[i]);
	err += fabs(gs[i]-fs[i]);
	size_d++;

	fs[i] = gs[i]/fs[i];
	}
	else
	{
	fs[i] = 1.0;
	}
	}
	info.volume = sum_f / sum_g;
	info.discrep = dis / size_d;
	info.error = err / size_d;*/

	return ErrorInfo();
}


void
Basis::
shrink(Matrix& c_err, const Matrix& c, const Matrix& l1, const Matrix& l2, fp shrinkage) const
{
	//cout << "  shrinkage Basis C_err" << c_err << " %*% C" << c << " %/% (" << shrinkage << " x L2" << l2 << " + L1" << l1 << " =";
	c_err.shrink(c, l1, l2, shrinkage);
	//cout << " C" << c_err << endl;

	// Standard shrinkage
	/*#pragma omp parallel for
	for (li i = 0; i < es.size(); i++)
	if (es[i] >= FLT_MIN && cs[i] >= FLT_MIN)
	{
	es[i] *= cs[i] / (shrinkage * l2[i] + wcs[i]);
	}
	else
	{
	es[i] = 0.0;
	}*/
}
