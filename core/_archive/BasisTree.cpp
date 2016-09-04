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


Basis::ErrorInfo
Basis::
error(SparseMatrix& fs, const SparseMatrix& gs) const
{
	ErrorInfo info;

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

	return info;
}


void
Basis::
shrink(SparseMatrix& es, const SparseMatrix& cs, const SparseMatrix& l2, const SparseMatrix& wcs, double shrinkage) const
{
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


////////////////////////////////////////////////////////////////////////////////


BasisBSpline::MeshInfo::
MeshInfo(ii _dc) :
	d(_dc),
	ls(_dc),
	os(_dc),
	ns(_dc),
	m(0)
{
}


void
BasisBSpline::MeshInfo::
operator=(const BasisBSpline::MeshInfo& ci)
{
	d = ci.d;
	ls = ci.ls;
	os = ci.os;
	ns = ci.ns;
	m = ci.m;
}


BasisBSpline::MeshInfo::
~MeshInfo()
{
}


li
BasisBSpline::MeshInfo::
size() const
{
	li size = m;
	for (ii i = 0; i < d; i++) size *= ns[i];
	return size;
}


void
BasisBSpline::MeshInfo::
print(ostream& out) const
{
	out << "l=[";
	for (ii i = 0; i < d; i++)
	{
		if (ls[i] == numeric_limits<ii>::min()) out << "?"; else out << ls[i];
		if (i < d - 1) out << ",";
	}
	out << "] oc=[";
	for (ii i = 0; i < d; i++)
	{
		if (os[i] == numeric_limits<ii>::min()) out << "?"; else out << os[i];
		if (i < d - 1) out << ",";
	}
	out << "] nc=" << m << "x[";
	for (ii i = 0; i < d; i++)
	{
		if (ns[i] == numeric_limits<ii>::min()) out << "?"; else out << ns[i];
		if (i < d - 1) out << ",";
	}
	out << "]=" << size();
}


