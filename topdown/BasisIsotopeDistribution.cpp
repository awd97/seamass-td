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


#include "BasisIsotopeDistribution.hpp"
#include <limits>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;


BasisIsotopeDistribution::
BasisIsotopeDistribution(vector<Basis*>& bases, Basis* parent, bool transient) :
	Basis(bases, parent->get_cm().d, parent, transient)
{
	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	cm = get_parent()->get_cm();

	m = parent->get_cm().n[0];

	nnz = 0;
	for (ii i = 0; i < m; i++)
	{
		if (i > 32) nnz++;
		nnz++;
	}
	vector<fp> acoo(nnz);
	vector<ii> rowind(nnz);
	vector<ii> colind(nnz);

	ii k = 0;
	for (ii i = 0; i < m; i++)
	{
		if (i > 32)
		{
			acoo[k] = 0.5;
			rowind[k] = i - 32;
			colind[k] = i;
			k++;
		}

		acoo[k] = 0.5;
		rowind[k] = i;
		colind[k] = i;
		k++;
	}

	a.resize(nnz);
	ia.resize(m + 1);
	ja.resize(nnz);

	ii job[] = { 2, 0, 0, 0, nnz, 0 }; ii info;
	mkl_scsrcoo(job, &m, a.data(), ja.data(), ia.data(), &nnz, acoo.data(), rowind.data(), colind.data(), &info);

	cout << index << " BasisIsotopeDistribution ";
	cm.print(cout);
	li size = sizeof(this);
	size += sizeof(fp) * a.capacity();
	size += sizeof(ii) * ia.capacity();
	size +=  sizeof(ii) * ja.capacity();
	cout << " mem=" << setprecision(2) << fixed << size / (1024.0*1024.0) << "Mb";
	if (transient) cout << " (t)";
	cout << endl;
	cout << "    A=[" << m << "," << cm.n[0] << "]:" << nnz << endl;
}


BasisIsotopeDistribution::~BasisIsotopeDistribution()
{
}


void
BasisIsotopeDistribution::
synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum)
{
	static fp alpha = 1.0;
	fp beta = accum ? 1.0 : 0.0;
	for (li j = 0; j < cm.n[1]; j++)
	{
		fp* c = const_cast<fp*>(&(cs.data()[j*cm.n[0]]));
		mkl_scsrmv("N", &m, &(cm.n[0]), &alpha, "G**C", a.data(), ja.data(), ia.data(), &(ia.data()[1]), c, &beta, &(fs.data()[j*cm.n[0]]));
	}
}


void
BasisIsotopeDistribution::
analysis(vector<fp>& es, const vector<fp>& fs)
{
	static fp alpha = 1.0, beta = 0.0;
	for (li j = 0; j < cm.n[1]; j++)
	{
		fp* f = const_cast<fp*>(&(fs.data()[j*cm.n[0]]));
		mkl_scsrmv("T", &m, &(cm.n[0]), &alpha, "G**C", a.data(), ja.data(), ia.data(), &(ia.data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
	}
}


void
BasisIsotopeDistribution::
l2norm(vector<fp>& es, const vector<fp>& fs)
{
	static fp alpha = 1.0, beta = 0.0;
	for (li j = 0; j < cm.n[1]; j++)
	{
		vector<fp> a2(nnz);
		vsSqr(nnz, a.data(), a2.data());
		fp* f = const_cast<fp*>(&(fs.data()[j*cm.n[0]]));
		mkl_scsrmv("T", &m, &(cm.n[0]), &alpha, "G**C", a2.data(), ja.data(), ia.data(), &(ia.data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
	}
}


void
BasisIsotopeDistribution::
shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage)
{
	parent->shrink(es, cs, l2, wcs, shrinkage);
}
