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


/*CoeffsMetadata::
CoeffsMetadata(ii _dc) :
	d(_dc), l(_dc), o(_dc), n(_dc)
{
}


void
CoeffsMetadata::
operator=(const CoeffsMetadata& cm)
{
	d = cm.d;
	l = cm.l;
	o = cm.o;
	n = cm.n;
}


CoeffsMetadata::
~CoeffsMetadata()
{
}


li
CoeffsMetadata::
size() const
{
	li size = 1;
	for (ii i = 0; i < d; i++) size *= n[i];
	return size;
}


void
CoeffsMetadata::
print(ostream& out) const
{
	out << "lc=[";
	for (ii i = 0; i < d; i++)
	{
		if (l[i] == numeric_limits<ii>::min()) out << "?"; else out << l[i];
		if (i < d - 1) out << ",";
	}
	out << "] oc=[";
	for (ii i = 0; i < d; i++)
	{
		if (o[i] == numeric_limits<ii>::min()) out << "?"; else out << o[i];
		if (i < d - 1) out << ",";
	}
	out << "] nc=[";
	for (ii i = 0; i < d; i++)
	{
		if (n[i] == numeric_limits<ii>::min()) out << "?"; else out << n[i];
		if (i < d - 1) out << ",";
	}
	out << "]:" << size();
}*/


Basis::
Basis(vector<Basis*>& bases, Basis* _parent, bool _transient) :
    parent(_parent),
    transient(_transient)
{
    index = bases.size();
    bases.push_back(this);
    if (parent) parent->child_count++;
}


Info
Basis::
error(vector<fp>& fs, const vector<fp>& gs) const
{
	Info info;

    ii size_d = 0;
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
    info.error = err / size_d;

	return info;
}


void
Basis::
shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage) const
{
	// Standard shrinkage
	#pragma omp parallel for
	for (li i = 0; i < es.size(); i++)
	if (es[i] >= FLT_MIN && cs[i] >= FLT_MIN)
	{
		es[i] *= cs[i] / (shrinkage * l2[i] + wcs[i]);
	}
	else
	{
		es[i] = 0.0;
	}
}
