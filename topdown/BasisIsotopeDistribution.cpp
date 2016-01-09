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
#include <map>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;


BasisIsotopeDistribution::
BasisIsotopeDistribution(vector<Basis*>& bases, Basis* parent,
                         ii out_res, ii factor_res, ii max_z,
                         bool transient) :
	Basis(bases, parent->get_cm().d, parent, transient)
{
	cout << index << " BasisIsotopeDistribution " << flush;

	ifstream ifs("seamass-td.factors");
	//ofstream ofs("grr");
	multimap< ii, vector<fp> > factors;
	ii i0 = parent->get_cm().o[0] / (1 << (out_res - factor_res));
	ii i1 = (parent->get_cm().o[0] + parent->get_cm().n[0] / max_z) / (1 << (out_res - factor_res));
	//cout << "i:[" << i0 << "," << i1 << "]" << endl;

	while (ifs.good())
	{
		char comma;
		ii i, z, id;
		ifs >> i;
		if (i < i0) { ifs.ignore(65535, '\n'); continue; }
		if (i > i1) break;
		
		ifs >> comma >> z >> comma >> id;
		map< ii, vector<fp> >::iterator elem = factors.insert(pair< ii, vector<fp> >(i, vector<fp>()));
		for (ii j = 0; j < 21; j++)
		{
			fp isotope;
			ifs >> comma >> isotope;
			if (elem->second.size() > 0 && isotope <= elem->second.back() && isotope < 0.00001)
			{
				ifs.ignore(65535, '\n');
				break;
			}
			else
			{
				elem->second.push_back(isotope);
			}
			//cout << " " << isotope << endl;
		}
		//ofs << i << " " << z << " " << id << " " << elem->second.size() << endl;

		if (i % 10000 == 0)
		{
			for (int i = 0; i < 256; ++i) cout << '\b';
			cout << index << " BasisIsotopeDistribution " << setw(1 + (int)(log10((float)i1))) << i0 << "/" << i << "/" << i1 << " " << flush;
		}
	}
	for (int i = 0; i < 256; ++i) cout << '\b';
	cout << index << " BasisIsotopeDistribution " << flush;


	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	cm = get_parent()->get_cm();

	m = parent->get_cm().n[0];

	nnz = m;
	vector<fp> acoo(nnz);
	vector<ii> rowind(nnz);
	vector<ii> colind(nnz);

	ii k = 0;
	for (ii j = 0; j < cm.n[1]; j++)
	for (ii i = 0; i < cm.n[0]; i++)
	{
		//ii f = (get_cm().o[0] + i) / (1 << (out_res - factor_res));
		//find->

		//for (ii z = 0; z < max_z; z++)
		{
			acoo[k] = 1.0;
			rowind[k] = k;
			colind[k] = k;
			k++;
		}
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
