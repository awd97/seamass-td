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


#include "BasisBSpline.hpp"
using namespace std;


BasisBSpline::
BasisBSpline(std::vector<Basis*>& bases, ii dims, bool transient, ii _parent_index) :
	Basis(bases, transient, _parent_index),
	mi(dims)
{
}

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


ii
BasisBSpline::MeshInfo::
n() const
{
	ii n = 1;
	for (ii i = 0; i < d; i++) n *= ns[i];
	return n;
}


li
BasisBSpline::MeshInfo::
size() const
{
	li size = m;
	for (ii i = 0; i < d; i++) size *= ns[i];
	return size;
}


ostream&
operator<<(ostream& os, const BasisBSpline::MeshInfo& mi)
{
	os << "[" << mi.m << ",";
	if (mi.d > 1) os << "(";
	for (ii i = 0; i < mi.d; i++)
	{
		os << mi.ns[i];
		if (i < mi.d - 1) os << ",";
	}
	if (mi.d > 1) os << ")";
	os << "]";

	/*out << "l=[";
	for (ii i = 0; i < mi.d; i++)
	{
		if (mi.ls[i] == numeric_limits<ii>::min()) out << "?"; else out << mi.ls[i];
		if (i < mi.d - 1) out << ",";
	}
	out << "] oc=[";
	for (ii i = 0; i < mi.d; i++)
	{
		if (mi.os[i] == numeric_limits<ii>::min()) out << "?"; else out << mi.os[i];
		if (i < mi.d - 1) out << ",";
	}
	out << "] nc=" << mi.m << "x[";
	for (ii i = 0; i < mi.d; i++)
	{
		if (mi.ns[i] == numeric_limits<ii>::min()) out << "?"; else out << mi.ns[i];
		if (i < mi.d - 1) out << ",";
	}
	out << "]=" << mi.size();*/

	return os;
}


